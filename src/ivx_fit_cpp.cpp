//' @export
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List ivx_fit_cpp(const arma::vec & y, const arma::mat & X, int K = 1) {

  int nr = X.n_rows;

  arma::mat xlag = X.rows(0, nr-2);
  arma::mat xt = X.rows(1, nr -1);
  arma::colvec yt = y.rows(1, nr -1);

  int  nn = xlag.n_rows, l = xlag.n_cols;
  // NumericVector df(2); df(0) = l, df(1) = nn -l;

  //join_horiz to include intercept
  arma::mat Xols = join_rows(ones(nr-1, 1), xlag);

  arma::colvec Aols = arma::solve(Xols, yt); //inv(x_con.t()*x_con)*x_con.t()*yt;
  arma::colvec epshat = yt - Xols*Aols;

  double s2 = std::inner_product(epshat.begin(), epshat.end(), epshat.begin(), 0.0)/(nn - l);
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(Xols)*Xols)));
  arma::colvec tstat = Aols/std_err;

  arma::mat Rn = zeros<mat>(l, l);
  for (int i = 0; i < l; i++) {
    double coef_Rn = as_scalar(arma::solve(xlag.col(i), xt.col(i)));
    Rn(i,i) = coef_Rn;
  }

  // // autoregressive residual estimation
  arma::mat u = xt - xlag * Rn;

  // // residuals' correlation matrix
  arma::mat corrmat = cor(epshat,u);

  // //covariance matrix estimation (predictive regression)
  arma::vec covepshat = epshat.t()*epshat/nn;

  // // covariance matrix estimation (autoregression)
  arma::mat covu = zeros<mat>(l, l);
  for (int i = 0; i < nn; ++i) {
    covu += trans(u.row(i))*u.row(i);
  }
  covu = covu / nn;

  // // covariance matrix between 'epshat' and 'u'
  arma::mat covuhat = zeros<mat>(1, l);
  for (int i = 0; i < l; ++i) {
    covuhat(i) = sum(epshat % u.col(i));
  }
  covuhat = covuhat.t() / nn;

  double m = floor(pow(nn, 0.3333333));
  arma::mat uu = zeros<mat>(l,l);
  for (int h = 1; h <= m; ++h) {
    arma::mat a = zeros<mat>(l,l);
    for (int t = h; t < nn; ++t) {
      a += u.row(t).t()*u.row(t-h);
    }
    // a.print("a:");
    double con = as_scalar(1 - h/(1+m));
    uu += con*a;
  }
  uu = uu/nn;
  arma::mat Omegauu = covu + uu + uu.t();

  arma::mat q = zeros<mat>(m,l);
  for (int h = 1; h <= m; ++h){
    arma::mat p = zeros<mat>(nn-h,l);
    for (int t = h; t < nn; ++t){
      p.row(t-h) = u.row(t) * as_scalar(epshat.row(t-h)); //1x1 matrix reduce to scalar
    }
    double con = as_scalar(1 - h/(1+m));
    q.row(h-1)= con*sum(p);
  }
  arma::mat residue = sum(q)/nn;
  arma::mat Omegaeu = covuhat + residue.t();

  // instrument construction
  arma::mat Rz = (1-1/(pow(nn, 0.95)))*eye(l,l);
  arma::mat diffx = xt - xlag;
  arma::mat z = zeros<mat>(nn,l);
  z.row(0) = diffx.row(0);
  for (int i = 1; i<nn; ++i){
    z.row(i) = z.row(i-1)*Rz + diffx.row(i);
  }

  int n = nn - K + 1;
  arma::mat Z = join_vert(zeros<mat>(1,l), z.rows(0,n-2));
  arma::mat zz = join_vert(zeros<mat>(1,l), z.rows(0,nn-2));

  arma::mat ZK = zeros<mat>(n,l);
  for (int i = 0; i < n; ++i){
    ZK.row(i) = sum(zz.rows(i, i+K-1)); // here should be sum(..., 1)
  }
  arma::mat meanzK = mean(ZK);

  arma::vec yy = zeros<vec>(n);
  for (int i = 0; i < n; ++i){
    yy.row(i)=sum(yt.rows(i, i+K-1));
  }
  arma::vec Yt = yy - as_scalar(mean(yy, 0));

  arma::mat xK = zeros<mat>(n,l);
  for (int i = 0; i<n; ++i){
    xK.row(i)=sum(xlag.rows(i,i+K-1));
  }
  arma::mat meanxK = mean(xK);

  arma::mat Xt = zeros<mat>(n,l);
  for (int i = 0; i < l; ++i){
    Xt.col(i) = xK.col(i) - ones(n,1)*meanxK.col(i);
  }

  ////////////////////////////////////////////////

  arma::mat Aivx = Yt.t()*Z*pinv(Xt.t()*Z);
  arma::mat intercept = mean(Yt) - mean(Xt) * Aivx.t();
  arma::colvec fitted = Xt*trans(Aivx);
  arma::colvec residuals = Yt - fitted;

  arma::mat FM = covepshat - Omegaeu.t()*inv(Omegauu)*Omegaeu;
  arma::mat M = ZK.t()*ZK*as_scalar(covepshat)-n*meanzK.t()*meanzK*as_scalar(FM);
  arma::mat H = eye<mat>(l,l);
  arma::mat Q = H*pinv(Z.t()*Xt)*M*pinv(Xt.t()*Z)*H.t();

  arma::colvec wivx = (H*Aivx.t()).t()*pinv(Q)*(H*Aivx.t());

  arma::mat wivxind_z = Aivx/sqrt(diagvec(Q).t());
  arma::mat wivxind = pow(wivxind_z.t(), 2);

  ////////////////////////////////////////////////

  return List::create(
    _("Aivx") = Aivx.t(),
    _("intercept") = intercept,
    _("fitted") = fitted,
    _("residuals") = residuals,
    _("wivx") = wivx,
    _("wivxind") = wivxind,
    _("zinvxind") = wivxind_z,
    _("Aols") = Aols,
    _("tstat_ols") = tstat,
    _("residuals_ols") = epshat,
    _("rank") = rank(Xt),
    _("rank_ols") = rank(X),
    _("horizons") = K,
    _("df.residuals") = nn - l,
    _("df") = l,
    _("delta") = corrmat,
    _("Rn") = diagvec(Rn),
    _("Rz") = diagvec(Rz),
    _("varcov") = Q
  );

}

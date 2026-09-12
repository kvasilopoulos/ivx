//' @export
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// rolling sum over K consecutive rows: row i = sum(A.rows(i, i+K-1))
static arma::mat roll_sum(const arma::mat & A, int K) {
  int n = A.n_rows - K + 1;
  arma::mat out(n, A.n_cols);
  for (int i = 0; i < n; ++i) out.row(i) = sum(A.rows(i, i+K-1), 0);
  return out;
}

// [[Rcpp::export]]
List ivx_fit_cpp(const arma::vec & y, const arma::mat & X, int K = 1,
                 double beta = 0.95, double cz = 1, int bandwidth = -1,
                 bool robust = false) {

  if (robust && K != 1) stop("robust = TRUE is only available for horizon = 1");

  const int nr = X.n_rows;
  const arma::mat xlag = X.rows(0, nr-2);
  const arma::mat xt = X.rows(1, nr-1);
  const arma::vec yt = y.rows(1, nr-1);
  const int nn = xlag.n_rows, l = xlag.n_cols;

  ///////////// OLS with intercept ///////////////

  arma::mat Xols = join_rows(ones(nn, 1), xlag);
  arma::vec Aols = solve(Xols, yt);
  arma::vec epshat = yt - Xols*Aols;
  double s2 = dot(epshat, epshat)/(nn - l);
  arma::vec std_err = sqrt(s2 * diagvec(pinv(Xols.t()*Xols)));
  arma::vec tstat = Aols/std_err;

  ///////// AR(1) per regressor, no intercept ////

  arma::vec Rn(l);
  for (int i = 0; i < l; ++i)
    Rn(i) = dot(xlag.col(i), xt.col(i)) / dot(xlag.col(i), xlag.col(i));
  arma::mat u = xt - xlag.each_row() % Rn.t();

  arma::mat corrmat = cor(epshat, u);
  double covepshat = dot(epshat, epshat)/nn;
  arma::mat covu = u.t()*u/nn;
  arma::vec covuhat = u.t()*epshat/nn;

  // Newey-West long-run covariances; default bandwidth n^(1/3) as in KMS (2015)
  int m = bandwidth < 0 ? (int) floor(pow(nn, 0.3333333)) : bandwidth;
  arma::mat uu(l, l, fill::zeros);
  arma::vec residue(l, fill::zeros);
  for (int h = 1; h <= std::min(m, nn-1); ++h) {
    double con = 1 - (double) h/(1+m);
    arma::mat ut = u.rows(h, nn-1).t();
    uu += con * (ut * u.rows(0, nn-h-1));
    residue += con * (ut * epshat.rows(0, nn-h-1));
  }
  arma::mat Omegauu = covu + (uu + uu.t())/nn;
  arma::vec Omegaeu = covuhat + residue/nn;

  ////////// instrument construction ////////////

  // instrument persistence: Rz = (1 - cz/n^beta) I, KMS use beta = 0.95, cz = 1
  double rz = 1 - cz/pow(nn, beta);
  arma::mat diffx = xt - xlag;
  arma::mat z(nn, l);
  z.row(0) = diffx.row(0);
  for (int i = 1; i < nn; ++i) z.row(i) = rz*z.row(i-1) + diffx.row(i);

  int n = nn - K + 1;
  arma::mat zz = join_vert(zeros<mat>(1, l), z.rows(0, nn-2)); // lagged instrument
  arma::mat Z = zz.rows(0, n-1);

  arma::mat ZK = roll_sum(zz, K);
  arma::rowvec meanzK = mean(ZK);

  arma::vec yy = roll_sum(yt, K);
  arma::vec Yt = yy - mean(yy);

  arma::mat xK = roll_sum(xlag, K);
  arma::mat Xt = xK.each_row() - mean(xK);

  ////////////////////////////////////////////////

  arma::mat XZinv = pinv(Xt.t()*Z);
  arma::rowvec Aivx = Yt.t()*Z*XZinv;
  arma::vec fitted = Xt*Aivx.t();
  arma::mat intercept = mean(Yt) - mean(Xt)*Aivx.t();
  arma::vec residuals = Yt - fitted;

  ///////////////// No demeaning /////////////////
  arma::mat interceptm = mean(y) - mean(xlag)*Aivx.t();
  arma::vec fittedm = as_scalar(interceptm) + xlag*Aivx.t();

  double FM = covepshat - as_scalar(Omegaeu.t()*solve(Omegauu, Omegaeu));
  // Eicker-White form (Demetrescu et al. 2023, Remarks 8-9): sigma^2 Z'Z -> sum z z' u^2
  arma::mat ZZ;
  if (robust) ZZ = ZK.t()*(ZK.each_col() % square(epshat));
  else ZZ = ZK.t()*ZK*covepshat;
  arma::mat M = ZZ - n*meanzK.t()*meanzK*FM;
  arma::mat Q = XZinv.t()*M*XZinv;

  arma::mat wivx = Aivx*pinv(Q)*Aivx.t();
  arma::mat wivxind_z = Aivx/sqrt(diagvec(Q).t());
  arma::mat wivxind = square(wivxind_z.t());

  ////////////////////////////////////////////////

  List data = List::create(
    _("X") = xlag,
    _("y") = yt.col(0)
    );

  List initial = List::create(
    _("intercept") = interceptm,
    _("fitted") = fittedm
    );

  List datam = List::create(
    _("X") = Xt,
    _("y") = Yt.col(0)
  );

  List ols = List::create(
    _("Aols") = Aols,
    _("se") = std_err,
    _("tstat_ols") = tstat,
    _("residuals_ols") = epshat,
    _("rank_ols") = rank(X)
    );

  return List::create(
    _("Aivx") = Aivx.t(),
    _("intercept") = intercept,
    _("fitted") = fitted,
    _("residuals") = residuals,
    _("wivx") = wivx,
    _("wivxind") = wivxind,
    _("zinvxind") = wivxind_z,
    _("rank") = rank(Xt),
    _("horizons") = K,
    _("df.residuals") = nn - l,
    _("df") = l,
    _("delta") = corrmat,
    _("Rn") = Rn,
    _("Rz") = arma::vec(l, fill::value(rz)),
    _("varcov") = Q,
    _("bandwidth") = m,
    _("ols") = ols,
    _("data") = data,
    _("initial") = initial,
    _("datam") = datam
  );

}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

static arma::mat roll_sum_sys(const arma::mat & A, int K) {
  int n = A.n_rows - K + 1;
  arma::mat out(n, A.n_cols);
  for (int i = 0; i < n; ++i) out.row(i) = sum(A.rows(i, i+K-1), 0);
  return out;
}

// Systems IVX (Kostakis, Magdalinos & Stamatogiannis 2023, eqs 15 and 23):
// Y is n x m (m responses), X is n x l (l predictors), A is m x l.
// Covariance of vec(A) is [(Z'X)^-1 (x) I_m] M [(X'Z)^-1 (x) I_m] with
// M = ZK'ZK (x) Sigma - n zbar zbar' (x) Sigma_FM. Reduces to ivx_fit_cpp for m = 1.
// [[Rcpp::export]]
List ivx_sys_fit_cpp(const arma::mat & Y, const arma::mat & X, int K = 1,
                     double beta = 0.95, double cz = 1, int bandwidth = -1) {
  const int nr = X.n_rows;
  const arma::mat xlag = X.rows(0, nr-2);
  const arma::mat xt = X.rows(1, nr-1);
  const arma::mat yt = Y.rows(1, nr-1);
  const int nn = xlag.n_rows, l = xlag.n_cols, m = yt.n_cols;

  // OLS with intercept, equation by equation
  arma::mat Xols = join_rows(ones(nn, 1), xlag);
  arma::mat Aols = solve(Xols, yt);
  arma::mat epshat = yt - Xols*Aols;                    // nn x m
  arma::mat Sigma = epshat.t()*epshat/nn;               // m x m

  // AR(1) per regressor, no intercept
  arma::vec Rn(l);
  for (int i = 0; i < l; ++i)
    Rn(i) = dot(xlag.col(i), xt.col(i)) / dot(xlag.col(i), xlag.col(i));
  arma::mat u = xt - xlag.each_row() % Rn.t();          // nn x l
  arma::mat corrmat = cor(epshat, u);                    // m x l

  // Newey-West long-run covariances
  int mb = bandwidth < 0 ? (int) floor(pow(nn, 0.3333333)) : bandwidth;
  arma::mat covu = u.t()*u/nn;
  arma::mat covue = u.t()*epshat/nn;                     // l x m
  arma::mat uu(l, l, fill::zeros), residue(l, m, fill::zeros);
  for (int h = 1; h <= std::min(mb, nn-1); ++h) {
    double con = 1 - (double) h/(1+mb);
    arma::mat ut = u.rows(h, nn-1).t();
    uu += con * (ut * u.rows(0, nn-h-1));
    residue += con * (ut * epshat.rows(0, nn-h-1));
  }
  arma::mat Omegauu = covu + (uu + uu.t())/nn;
  arma::mat Omegaue = covue + residue/nn;                // l x m
  arma::mat FM = Sigma - Omegaue.t()*solve(Omegauu, Omegaue);   // m x m

  // instrument
  double rz = 1 - cz/pow(nn, beta);
  arma::mat diffx = xt - xlag;
  arma::mat z(nn, l);
  z.row(0) = diffx.row(0);
  for (int i = 1; i < nn; ++i) z.row(i) = rz*z.row(i-1) + diffx.row(i);

  int n = nn - K + 1;
  arma::mat zz = join_vert(zeros<mat>(1, l), z.rows(0, nn-2));
  arma::mat Z = zz.rows(0, n-1);
  arma::mat ZK = roll_sum_sys(zz, K);
  arma::rowvec meanzK = mean(ZK);

  arma::mat yy = roll_sum_sys(yt, K);
  arma::mat Yt = yy.each_row() - mean(yy);
  arma::mat xK = roll_sum_sys(xlag, K);
  arma::mat Xt = xK.each_row() - mean(xK);

  // estimator (15): A = Y'Z (X'Z)^-1, m x l
  arma::mat XZinv = pinv(Xt.t()*Z);
  arma::mat A = Yt.t()*Z*XZinv;
  arma::mat fitted = Xt*A.t();
  arma::mat residuals = Yt - fitted;
  arma::rowvec intercept = mean(Yt) - mean(Xt)*A.t();

  // covariance (23) of vec(A)
  arma::mat Im = eye(m, m);
  arma::mat M = kron(ZK.t()*ZK, Sigma) - n*kron(meanzK.t()*meanzK, FM);
  arma::mat G = kron(XZinv.t(), Im);                    // (Z'X)^-1 (x) I_m
  arma::mat Q = G*M*G.t();
  arma::vec a = vectorise(A);
  double wald = as_scalar(a.t()*solve(Q, a));
  arma::mat se = reshape(sqrt(Q.diag()), m, l);
  arma::mat tstat = A/se;

  // per-equation Wald: H selects the m-spaced entries of vec(A)
  arma::vec wald_eq(m);
  for (int i = 0; i < m; ++i) {
    arma::uvec idx = regspace<uvec>(i, m, m*l - 1);
    arma::vec ai = a.elem(idx);
    arma::mat Qi = Q.submat(idx, idx);
    wald_eq(i) = as_scalar(ai.t()*solve(Qi, ai));
  }

  return List::create(
    _("A") = A,
    _("se") = se,
    _("tstat") = tstat,
    _("intercept") = intercept,
    _("fitted") = fitted,
    _("residuals") = residuals,
    _("wald") = wald,
    _("wald_eq") = wald_eq,
    _("df") = m*l,
    _("df.residuals") = nn - l,
    _("delta") = corrmat,
    _("Rn") = Rn,
    _("Rz") = rz,
    _("varcov") = Q,
    _("bandwidth") = mb,
    _("ols") = List::create(_("Aols") = Aols, _("residuals") = epshat, _("Sigma") = Sigma)
  );
}

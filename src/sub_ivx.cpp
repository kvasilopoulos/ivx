// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Subsample IVX statistics (Demetrescu et al. 2023, eqs 15-17 and Remark 11).
// y: y_t, x: x_{t-1}, z: z_{t-1} (full-sample instrument), rows aligned;
// windows are [start_i, end_i] (0-based, inclusive). Returns one row per
// window: Wald statistic followed by the t-ratio of each coefficient.
// [[Rcpp::export]]
arma::mat sub_ivx_cpp(const arma::vec & y, const arma::mat & x, const arma::mat & z,
                      const arma::uvec & start, const arma::uvec & end, bool robust) {
  int l = x.n_cols, nw = start.n_elem;
  arma::mat out(nw, 1 + l);
  for (int i = 0; i < nw; ++i) {
    arma::vec ys = y.rows(start(i), end(i));
    arma::mat xs = x.rows(start(i), end(i));
    arma::mat zs = z.rows(start(i), end(i));
    ys -= arma::mean(ys);
    xs.each_row() -= arma::mean(xs, 0);
    arma::mat B = zs.t() * xs;
    arma::vec b = arma::solve(B, zs.t() * ys);
    arma::vec u = ys - xs * arma::solve(xs, ys);          // subsample OLS residuals
    arma::mat M;
    if (robust) M = zs.t() * arma::diagmat(arma::square(u)) * zs;
    else M = zs.t() * zs * arma::mean(arma::square(u));
    arma::mat Binv = arma::inv(B);
    arma::mat V = Binv * M * Binv.t();
    out(i, 0) = arma::as_scalar(b.t() * arma::solve(V, b));
    for (int k = 0; k < l; ++k) out(i, 1 + k) = b(k) / std::sqrt(V(k, k));
  }
  return out;
}

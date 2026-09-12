// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// x_t = sum_j A_j x_{t-j} + v_t with zero initial conditions
// (Demetrescu et al. 2023, Algorithm 1, Step 4). A is a list of l x l lag matrices.
// [[Rcpp::export]]
arma::mat var_sim_cpp(const arma::mat & v, const List & A) {
  int n = v.n_rows, q = A.size();
  std::vector<arma::mat> As(q);
  for (int j = 0; j < q; ++j) As[j] = as<arma::mat>(A[j]);
  arma::mat x = v;
  for (int t = 1; t < n; ++t) {
    for (int j = 0; j < std::min(q, t); ++j) {
      x.row(t) += x.row(t - 1 - j) * As[j].t();
    }
  }
  return x;
}

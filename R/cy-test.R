#' Bonferroni Q-Test of Campbell and Yogo (2006)
#'
#' `cy_test` implements the Bonferroni Q-test of Campbell and Yogo (2006), the
#' standard feasible version of the sup-bound / Bonferroni approach of
#' Cavanagh, Elliott and Stock (1995). For a single predictor
#' \eqn{x_t = \gamma + \rho x_{t-1} + v_t} (AR(p) dynamics allowed, Appendix A),
#' the procedure is:
#' \enumerate{
#' \item estimate the innovation correlation \eqn{\delta} between the
#'   predictive-regression residual and the ADF innovation of the predictor
#'   (lag length by BIC, \eqn{p \in [1, p_{max}]});
#' \item compute the DF-GLS statistic of Elliott, Rothenberg and Stock (1996)
#'   and invert its local-to-unity null distribution (Stock, 1991) into a
#'   confidence interval \eqn{[\underline c, \bar c]} for \eqn{c = T(\rho - 1)},
#'   with lower and upper levels \eqn{\underline\alpha_1(\delta)},
#'   \eqn{\bar\alpha_1(\delta)} from the paper's Table 2, which refine the
#'   Bonferroni bound so that the one-sided test has size 5\%;
#' \item for \eqn{\rho} at each end of the interval compute the Q-estimate
#'   \eqn{\hat\beta(\rho)} of eq. (25) - OLS of \eqn{y_t} on the demeaned
#'   \eqn{x_{t-1}} after subtracting \eqn{(\sigma_{ue}/\sigma_e\omega)(x_t - \rho x_{t-1})},
#'   with the Phillips-Perron-type correction \eqn{\tfrac{T}{2}(\sigma_{ue}/\sigma_e\omega)(\omega^2 - \sigma_v^2)}
#'   when \eqn{p > 1} - and its standard error \eqn{\sigma_u (1-\delta^2)^{1/2}/(\sum x^{\mu 2}_{t-1})^{1/2}};
#' \item the 90\% Bonferroni confidence interval for \eqn{\beta} is
#'   \eqn{[\hat\beta(\bar\rho) - 1.645\,se,\; \hat\beta(\underline\rho) + 1.645\,se]}
#'   (eq. 17); the 5\% one-sided Q-tests reject \eqn{\beta \le 0} if the lower
#'   bound is positive and \eqn{\beta \ge 0} if the upper bound is negative.
#' }
#' The DF-GLS null quantiles are tabulated by simulation for
#' \eqn{c \in [-100, 10]} (see `data-raw/dfgls-quantiles.R`); a statistic outside
#' the tabulated range is clamped to the boundary, which for very negative
#' values (a clearly stationary predictor) makes the interval for \eqn{c}
#' start at \eqn{-100}. Table 2 is given for \eqn{\delta < 0}; for
#' \eqn{\hat\delta > 0} the predictor is sign-flipped, which flips \eqn{\beta}
#' and the alternative, and the results are mapped back.
#'
#' @inheritParams ivx
#' @param lag_max maximum ADF lag order for the BIC search; the default is
#' \eqn{\lfloor 12 (T/100)^{1/4} \rfloor} lagged differences.
#'
#' @return an object of class "cy_test": a list with `ci` (the 90\% Bonferroni
#' interval for \eqn{\beta}), `reject` (named logical: `greater`, `less`),
#' `estimate` (OLS slope), `beta_rho` (\eqn{\hat\beta} at \eqn{\underline\rho}
#' and \eqn{\bar\rho}), `se`, `delta`, `dfgls`, `c_ci`, `rho_ci`, `alpha1`,
#' `lag` and `n`.
#'
#' @references Campbell, J. Y., & Yogo, M. (2006). Efficient tests of stock
#' return predictability. Journal of Financial Economics, 81(1), 27-60.
#' @references Cavanagh, C. L., Elliott, G., & Stock, J. H. (1995). Inference in
#' models with nearly integrated regressors. Econometric Theory, 11(5), 1131-1147.
#' @references Elliott, G., Rothenberg, T. J., & Stock, J. H. (1996). Efficient
#' tests for an autoregressive unit root. Econometrica, 64(4), 813-836.
#' @references Stock, J. H. (1991). Confidence intervals for the largest
#' autoregressive root in U.S. macroeconomic time series. Journal of Monetary
#' Economics, 28(3), 435-459.
#'
#' @export
#' @examples
#' cy_test(Ret ~ DP, data = kms)
cy_test <- function(formula, data, lag_max = NULL, na.action) {
  cl <- match.call()
  fr <- ivx_frame(match.call(expand.dots = FALSE), parent.frame())
  if (NCOL(fr$x) != 1) stop("cy_test() is defined for a single predictor", call. = FALSE)
  out <- cy_test_fit(fr$y, drop(fr$x), lag_max = lag_max)
  out$call <- cl
  out$predictor <- colnames(fr$x)
  out
}

#' @rdname cy_test
#' @param y response vector.
#' @param x predictor vector.
#' @export
cy_test_fit <- function(y, x, lag_max = NULL) {
  n <- length(x)
  if (length(y) != n) stop("incompatible dimensions")
  pmax <- if (is.null(lag_max)) floor(12 * (n / 100)^0.25) else as.integer(lag_max)

  # ADF regression of the predictor, k = p - 1 lagged differences chosen by BIC
  adf <- adf_lag(x, pmax, ic = "bic")
  p <- adf$k + 1
  e <- adf$resid                                    # innovations e_t, t = p + 1..n
  psi <- adf$psi
  # predictive regression on the same sample
  idx <- (p + 1):n
  xl <- x[idx - 1]
  yt <- y[idx]
  ols <- lm.fit(cbind(1, xl), yt)
  u <- ols$residuals
  Tn <- length(idx)
  s2u <- sum(u^2) / (Tn - 2)
  s2e <- sum(e^2) / (Tn - p - 1)   # p + 1 ADF coefficients
  sue <- sum(u * e) / Tn
  delta <- sue / sqrt(s2u * s2e)
  # long-run scale omega = sigma_e / b(1) and var(v_t), v_t the AR(p - 1) error
  b1 <- 1 - sum(psi)
  omega2 <- s2e / b1^2
  s2v <- if (p > 1) ar_var(psi, s2e) else s2e

  flip <- delta > 0
  if (flip) {
    x <- -x; xl <- -xl; e <- -e; sue <- -sue; delta <- -delta
  }

  # DF-GLS with the same lag order and the CI for c from CY Table 2 levels
  dfgls <- dfgls_stat(x, p)
  a1 <- cy_alpha1(delta)
  c_lo <- dfgls_invert(dfgls, 1 - a1[["a1_lower"]])
  c_hi <- dfgls_invert(dfgls, a1[["a1_upper"]])
  rho_lo <- 1 + c_lo / n
  rho_hi <- 1 + c_hi / n

  # Q-estimate of beta given rho (eq. 25) and its standard error (eqs 15-16)
  xm <- xl - mean(xl)
  sxx <- sum(xm^2)
  k <- sue / sqrt(s2e * omega2)
  beta_rho <- function(rho) {
    (sum(xm * (yt - k * (x[idx] - rho * xl))) - Tn / 2 * k * (omega2 - s2v)) / sxx
  }
  se <- sqrt(s2u * (1 - delta^2) / sxx)
  b_lo <- beta_rho(rho_hi) - stats::qnorm(0.95) * se   # lower bound at the upper rho
  b_hi <- beta_rho(rho_lo) + stats::qnorm(0.95) * se
  b_at <- c(rho_lower = beta_rho(rho_lo), rho_upper = beta_rho(rho_hi))
  ci <- c(lower = b_lo, upper = b_hi)
  est <- unname(ols$coefficients[2])
  if (flip) {
    ci <- c(lower = -b_hi, upper = -b_lo)
    b_at <- -b_at[2:1]
    names(b_at) <- c("rho_lower", "rho_upper")
    delta <- -delta
  }
  structure(
    list(
      ci = ci, reject = c(greater = unname(ci["lower"] > 0), less = unname(ci["upper"] < 0)),
      estimate = est, beta_rho = b_at, se = se, delta = delta, dfgls = dfgls,
      c_ci = c(lower = c_lo, upper = c_hi), rho_ci = c(lower = rho_lo, upper = rho_hi),
      alpha1 = a1, lag = p, n = n, flipped = flip
    ),
    class = "cy_test"
  )
}

# DF-GLS statistic (ERS 1996): quasi-GLS demeaning with c_bar = -7, then the ADF
# t-ratio on the demeaned series with p - 1 lagged differences and no intercept
dfgls_stat <- function(x, p) {
  n <- length(x)
  rb <- 1 - 7 / n
  z <- c(x[1], x[-1] - rb * x[-n])
  d <- c(1, rep(1 - rb, n - 1))
  xd <- x - sum(d * z) / sum(d^2)
  dx <- diff(xd)
  k <- p - 1
  idx <- (k + 2):n
  X <- cbind(xd[idx - 1])
  if (k > 0) X <- cbind(X, sapply(seq_len(k), function(i) dx[idx - 1 - i]))
  f <- lm.fit(X, dx[idx - 1])
  s2 <- sum(f$residuals^2) / (length(idx) - ncol(X))
  unname(f$coefficients[1] / sqrt(s2 * solve(crossprod(X))[1, 1]))
}

# c such that the p-quantile of the DF-GLS null distribution at c equals the
# statistic; linear interpolation on the simulated grid, clamped to its range
dfgls_invert <- function(stat, prob) {
  ps <- as.numeric(rownames(dfgls_q))
  cs <- as.numeric(colnames(dfgls_q))
  q <- sapply(seq_along(cs), function(j) stats::approx(ps, dfgls_q[, j], xout = prob, rule = 2)$y)
  q <- cummax(q)   # simulation noise; the quantile is increasing in c
  if (stat <= q[1]) return(cs[1])
  if (stat >= q[length(q)]) return(cs[length(cs)])
  keep <- !duplicated(q)
  stats::approx(q[keep], cs[keep], xout = stat)$y
}

# Campbell & Yogo (2006) Table 2 levels, interpolated in delta (delta < 0)
cy_alpha1 <- function(delta) {
  d <- pmin(pmax(delta, -0.999), -0.025)
  c(a1_lower = stats::approx(cy_table2[, "delta"], cy_table2[, "a1_lower"], xout = d)$y,
    a1_upper = stats::approx(cy_table2[, "delta"], cy_table2[, "a1_upper"], xout = d)$y)
}

# variance of an AR(k) process with coefficients psi and innovation variance s2
ar_var <- function(psi, s2) {
  k <- length(psi)
  A <- if (k == 1) matrix(psi) else rbind(psi, cbind(diag(k - 1), 0))
  S <- matrix(0, k, k); S[1, 1] <- s2
  G <- matrix(solve(diag(k^2) - kronecker(A, A), as.vector(S)), k, k)
  G[1, 1]
}

#' @rdname cy_test
#' @param x an object of class "cy_test".
#' @param digits minimal number of significant digits.
#' @param ... unused.
#' @export
print.cy_test <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat("\nBonferroni Q-test (Campbell & Yogo, 2006)\n\n")
  cat("delta = ", formatC(x$delta, digits = digits), ", DF-GLS = ", formatC(x$dfgls, digits = digits),
      " (p = ", x$lag, "), CI for c at levels (", formatC(x$alpha1[1], digits = 2), ", ",
      formatC(x$alpha1[2], digits = 2), "): [", formatC(x$c_ci[1], digits = digits), ", ",
      formatC(x$c_ci[2], digits = digits), "], rho: [", formatC(x$rho_ci[1], digits = 4), ", ",
      formatC(x$rho_ci[2], digits = 4), "]\n", sep = "")
  cat("OLS slope = ", formatC(x$estimate, digits = digits), "; Q-estimates at the ends of the rho interval: ",
      formatC(x$beta_rho[1], digits = digits), ", ", formatC(x$beta_rho[2], digits = digits), "\n", sep = "")
  cat("90% Bonferroni confidence interval for beta: [", formatC(x$ci[1], digits = digits), ", ",
      formatC(x$ci[2], digits = digits), "]\n", sep = "")
  cat("5% one-sided Q-tests: H1 beta > 0 ", if (x$reject["greater"]) "reject" else "do not reject",
      " H0; H1 beta < 0 ", if (x$reject["less"]) "reject" else "do not reject", " H0\n\n", sep = "")
  invisible(x)
}

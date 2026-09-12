#' Hybrid t-Test for Return Predictability (Harvey, Leybourne & Taylor)
#'
#' `hlt_test` implements the double-switching hybrid procedure \eqn{T_{hyb}} of
#' Harvey, Leybourne and Taylor (2021, Section 3.3) for a single predictor.
#' Two regression t-ratios are used: the standard OLS t-ratio \eqn{T}
#' (eq. 5) and the variant \eqn{\tilde T} in which the predictor is quasi-GLS
#' demeaned with \eqn{\bar c = 7} (eq. 7). Under weak persistence the standard
#' t-ratio is compared with normal critical values; under strong persistence
#' the limiting null distributions depend on the local-to-unity parameter and
#' on the innovation correlation \eqn{\rho_{xy}}, and the tests are run with
#' the paper's asymptotically conservative critical values (maximised over
#' \eqn{c}) obtained from the response surfaces in its Table 1.
#'
#' The procedure is: (1) if the ADF normalised-bias statistic
#' \eqn{T\hat\rho/(1 - \sum_i \hat\psi_i)} from an ADF regression with lag
#' length chosen by the MBIC of Ng and Perron (2001) is below
#' \eqn{-4\sqrt{T}}, the predictor is treated as weakly persistent and the
#' standard test \eqn{T_N} is used; (2) otherwise, for an upper-tail test,
#' \eqn{T} with critical value \eqn{cv(\hat\rho_{xy})} if
#' \eqn{\hat\rho_{xy} > -0.1} and \eqn{\tilde T} with \eqn{\tilde{cv}(\hat\rho_{xy})}
#' if \eqn{\hat\rho_{xy} < -0.1} (mirrored for lower-tail tests), where
#' \eqn{\hat\rho_{xy}} is the correlation of the ADF residuals with the
#' predictive-regression residuals.
#'
#' @inheritParams ivx
#' @param alternative direction of the one-sided test.
#' @param level significance level; one of 0.1, 0.05, 0.025, 0.01 (the
#' levels for which response surfaces are available).
#' @param lag_max maximum ADF lag order; the default is the paper's
#' \eqn{\lfloor 12 (T/100)^{1/4} \rfloor}.
#'
#' @return an object of class "hlt_test": a list with the selected `test`
#' (`"T_N"`, `"T_con"` or `"T~_con"`), its `statistic`, `cv` and `reject`
#' indicator, plus `t`, `t_gls` (both t-ratios), `adf`, `adf_lag`,
#' `rho_xy` and `estimate` (the OLS slope).
#'
#' @references Harvey, D. I., Leybourne, S. J., & Taylor, A. M. R. (2021).
#' Simple tests for stock return predictability with good size and power
#' properties. Journal of Econometrics, 224(1), 198-214.
#' @references Ng, S., & Perron, P. (2001). Lag length selection and the
#' construction of unit root tests with good size and power. Econometrica,
#' 69(6), 1519-1554.
#'
#' @export
#' @examples
#' hlt_test(Ret ~ DP, data = kms)
#' hlt_test(Ret ~ TBL, data = kms, alternative = "less", level = 0.1)
hlt_test <- function(formula, data, alternative = c("greater", "less"), level = 0.05,
                     lag_max = NULL, na.action) {
  alternative <- match.arg(alternative)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  if (NCOL(x) != 1) stop("hlt_test() is defined for a single predictor", call. = FALSE)
  out <- hlt_test_fit(y, drop(x), alternative = alternative, level = level, lag_max = lag_max)
  out$call <- cl
  out$predictor <- colnames(x)
  out
}

#' @rdname hlt_test
#' @param y response vector.
#' @param x predictor vector.
#' @export
hlt_test_fit <- function(y, x, alternative = "greater", level = 0.05, lag_max = NULL) {
  if (!level %in% c(0.1, 0.05, 0.025, 0.01)) {
    stop("`level` must be one of 0.1, 0.05, 0.025, 0.01", call. = FALSE)
  }
  n <- length(x)
  if (length(y) != n) stop("incompatible dimensions")
  xl <- x[-n]
  yt <- y[-1]
  yd <- yt - mean(yt)

  # T: OLS demeaning of both series (eq. 5)
  f1 <- lm.fit(cbind(1, xl), yt)
  t_ols <- unname(f1$coefficients[2] / sqrt(sum(f1$residuals^2) / (n - 3) / sum((xl - mean(xl))^2)))

  # T~: quasi-GLS demeaned predictor, c_bar = 7 (Elliott et al., 1996), eq. (7)
  rb <- 1 - 7 / n
  xq <- c(x[1], x[-1] - rb * x[-n])
  dq <- c(1, rep(1 - rb, n - 1))
  mu_x <- sum(dq * xq) / sum(dq^2)
  xg <- xl - mu_x
  f2 <- lm.fit(cbind(xg), yd)
  t_gls <- unname(f2$coefficients[1] / sqrt(sum(f2$residuals^2) / (n - 2) / sum(xg^2)))

  # ADF regression with MBIC lag selection; normalised-bias statistic
  pmax <- if (is.null(lag_max)) floor(12 * (n / 100)^0.25) else as.integer(lag_max)
  adf <- adf_mbic(x, pmax)
  rho_xy <- cor(adf$resid, f1$residuals[(length(f1$residuals) - length(adf$resid) + 1):length(f1$residuals)])

  weak <- adf$stat < -4 * sqrt(n)
  upper <- alternative == "greater"
  r <- if (upper) rho_xy else -rho_xy
  if (weak) {
    test <- "T_N"; stat <- t_ols; cv <- stats::qnorm(1 - level)
  } else if (r > -0.1) {
    test <- "T_con"; stat <- t_ols; cv <- hlt_cv(r, level, gls = FALSE)
  } else {
    test <- "T~_con"; stat <- t_gls; cv <- hlt_cv(r, level, gls = TRUE)
  }
  if (!upper) cv <- -cv
  reject <- if (upper) stat > cv else stat < cv

  structure(
    list(
      test = test, statistic = stat, cv = cv, reject = reject, level = level,
      alternative = alternative, t = t_ols, t_gls = t_gls, adf = adf$stat,
      adf_lag = adf$p, rho_xy = rho_xy, estimate = unname(f1$coefficients[2]), n = n
    ),
    class = "hlt_test"
  )
}

# Table 1 response surfaces: cv = sum_k a_k rho^k, k = 0..8
hlt_rs <- list(
  ols = cbind(
    "0.1"   = c(1.346, -0.819, 1.928, -0.402, -5.008, 0.825, 7.040, -0.470, -3.607),
    "0.05"  = c(1.707, -0.802, 2.314, -0.377, -6.970, 1.013, 10.279, -0.705, -5.417),
    "0.025" = c(2.004, -0.765, 1.947, -0.602, -5.131, 0.965, 6.692, -0.350, -3.154),
    "0.01"  = c(2.434, -0.726, 1.257, -0.736, -2.385, 1.448, 1.972, -0.762, -0.479)
  ),
  gls = cbind(
    "0.1"   = c(1.293, -0.242, -0.055, -0.316, 0.493, 0.401, -0.808, -0.200, 0.459),
    "0.05"  = c(1.648, -0.225, 0.323, -0.275, -1.447, 0.432, 2.603, -0.290, -1.581),
    "0.025" = c(1.950, -0.285, 0.200, -0.186, -0.559, -0.005, 0.224, 0.219, 0.236),
    "0.01"  = c(2.377, -0.382, -0.171, 0.414, 0.209, -0.984, -0.504, 0.693, 0.434)
  )
)

hlt_cv <- function(rho, level, gls) {
  a <- hlt_rs[[if (gls) "gls" else "ols"]][, as.character(level)]
  sum(a * rho^(0:8))
}

# ADF regression dx_t = mu + rho x_{t-1} + sum psi_i dx_{t-i} + e_t, lag by the MBIC of
# Ng & Perron (2001) on the common sample t = pmax + 2..n, OLS-demeaned data
# (Perron & Qu, 2007). Returns the normalised bias T rho / (1 - sum psi), p, residuals.
adf_mbic <- function(x, pmax) {
  n <- length(x)
  dx <- diff(x)
  fit <- function(p, common) {
    start <- if (common) pmax + 1 else p + 1          # index into dx
    idx <- start:(n - 1)
    X <- cbind(1, x[idx])
    if (p > 0) X <- cbind(X, sapply(seq_len(p), function(i) dx[idx - i]))
    f <- lm.fit(X, dx[idx])
    list(coef = f$coefficients, resid = f$residuals, xl = x[idx])
  }
  ic <- sapply(0:pmax, function(p) {
    f <- fit(p, common = TRUE)
    Te <- length(f$resid)
    s2 <- sum(f$resid^2) / Te
    tau <- f$coef[2]^2 * sum(f$xl^2) / s2
    log(s2) + log(Te) * (tau + p) / Te
  })
  p <- which.min(ic) - 1
  f <- fit(p, common = FALSE)
  psi <- if (p > 0) f$coef[-(1:2)] else 0
  list(stat = unname(n * f$coef[2] / (1 - sum(psi))), p = p, resid = f$resid)
}

#' @rdname hlt_test
#' @param x an object of class "hlt_test".
#' @param digits minimal number of significant digits.
#' @param ... unused.
#' @export
print.hlt_test <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat("\nHybrid predictability test (Harvey, Leybourne & Taylor, 2021)\n\n")
  cat("Selected test: ", x$test,
      switch(x$test, "T_N" = " (weak persistence, normal critical value)",
             "T_con" = " (OLS-demeaned t, conservative critical value)",
             "T~_con" = " (quasi-GLS-demeaned t, conservative critical value)"), "\n", sep = "")
  cat("statistic = ", formatC(x$statistic, digits = digits), ", ", x$level * 100, "% critical value = ",
      formatC(x$cv, digits = digits), " (alternative: beta ", if (x$alternative == "greater") ">" else "<",
      " 0): ", if (x$reject) "reject" else "do not reject", "\n", sep = "")
  cat("slope = ", formatC(x$estimate, digits = digits), ", t = ", formatC(x$t, digits = digits),
      ", t (quasi-GLS) = ", formatC(x$t_gls, digits = digits), ", ADF = ", formatC(x$adf, digits = digits),
      " (p = ", x$adf_lag, ", cutoff ", formatC(-4 * sqrt(x$n), digits = digits), "), rho_xy = ",
      formatC(x$rho_xy, digits = digits), "\n\n", sep = "")
  invisible(x)
}

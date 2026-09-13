#' Control-Function Predictability Test of Elliott (2011)
#'
#' `elliott_cf` runs the augmented predictive regression of Elliott (2011,
#' eq. 9),
#' \deqn{y_t = \alpha + \beta' x_{t-1} + \gamma' Z_t + \tilde u_t,\qquad
#' Z_t = (z_t', z_{t-1}', \dots, z_{t-q}')',}
#' where \eqn{z_t} are user-supplied stationary covariates that are
#' contemporaneously correlated with the innovations of the persistent
#' predictors \eqn{x_t} and of the response. If they absorb that correlation
#' (the "orthogonalising" condition, Section 5 of the paper), the Wald test of
#' \eqn{\beta = 0} has a standard chi-square limit whatever the persistence of
#' \eqn{x_t} (Theorem 2); without them it has the non-standard Elliott-Stock
#' (1994) distribution (Theorem 1). Unlike [arm()] and [ivx_ra()], which build
#' the control variable from the data, the covariates here come from the user
#' - the paper's example is predicting returns with the dividend-price ratio
#' using contemporaneous price-related variables.
#'
#' The remaining innovation correlation after augmentation is returned as a
#' diagnostic: it is the correlation between the regression residuals and the
#' residuals of an AR(1) of each predictor on the same covariates, and should
#' be close to zero for the test to be reliable.
#'
#' @inheritParams ivx
#' @param covariates a one-sided formula giving the orthogonalising
#' covariates \eqn{z_t} (contemporaneous with the response).
#' @param lags number of lags \eqn{q} of the covariates to include.
#' @param robust logical; if `TRUE` (default) Eicker-White standard errors are
#' used.
#'
#' @return an object of class "elliott_cf": a list with `coefficients` (on the
#' predictors), `se`, `tstat`, `Wald`, `df`, `p.value`, `gamma` (coefficients on
#' the covariates and their lags), `rho_resid` (remaining innovation
#' correlation per predictor) and the underlying `lm` fit.
#'
#' @references Elliott, G. (2011). A control function approach for testing the
#' usefulness of trending variables in predictive regressions and econometric
#' models. Journal of Econometrics, 164(1), 79-91.
#' @references Elliott, G., & Stock, J. H. (1994). Inference in time series
#' regression when the order of integration of a regressor is unknown.
#' Econometric Theory, 10(3-4), 672-700.
#'
#' @export
#' @examples
#' # the T-bill rate as covariate for the dividend-price ratio (illustration only)
#' elliott_cf(Ret ~ DP, ~ TBL, data = kms)
elliott_cf <- function(formula, covariates, data, lags = 0, robust = TRUE, na.action) {
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
  zf <- stats::update(covariates, ~ . - 1)
  z <- model.matrix(zf, if (missing(data)) environment(covariates) else data)
  if (NROW(z) != NROW(x)) stop("covariates and predictors must have the same length", call. = FALSE)
  out <- elliott_cf_fit(y, x, z, lags = lags, robust = robust)
  out$call <- cl
  out
}

#' @rdname elliott_cf
#' @param y response vector.
#' @param x matrix of predictors.
#' @param z matrix of covariates.
#' @export
elliott_cf_fit <- function(y, x, z, lags = 0, robust = TRUE) {
  x <- as.matrix(x)
  z <- as.matrix(z)
  n <- NROW(x)
  k <- NCOL(x)
  if (NROW(y) != n || NROW(z) != n) stop("incompatible dimensions")
  q <- as.integer(lags)
  if (q < 0) stop("`lags` must be non-negative", call. = FALSE)
  cn <- colnames(x)
  if (is.null(cn)) cn <- paste0("x", seq_len(k))
  zn <- colnames(z)
  if (is.null(zn)) zn <- paste0("z", seq_len(NCOL(z)))

  # Z_t = (z_t, z_{t-1}, ..., z_{t-q}), aligned with y_t and x_{t-1}, t = q + 2..n
  idx <- (q + 2):n
  Z <- do.call(cbind, lapply(0:q, function(j) z[idx - j, , drop = FALSE]))
  colnames(Z) <- as.vector(outer(zn, 0:q, function(a, j) ifelse(j == 0, a, paste0(a, "_lag", j))))
  X <- cbind(x[idx - 1, , drop = FALSE], Z)
  colnames(X)[seq_len(k)] <- cn
  fit <- stats::lm(y[idx] ~ X)
  cf <- stats::coef(fit)[-1]
  names(cf) <- colnames(X)
  V <- if (robust) hc0(fit) else stats::vcov(fit)
  V <- V[-1, -1, drop = FALSE]
  beta <- cf[seq_len(k)]
  Vb <- V[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(Vb) <- list(cn, cn)
  se <- sqrt(diag(Vb))
  W <- drop(crossprod(beta, solve(Vb, beta)))

  # remaining innovation correlation: residuals of x_t on x_{t-1} and Z_t
  e <- stats::residuals(fit)
  rho <- sapply(seq_len(k), function(j) {
    ex <- lm.fit(cbind(1, x[idx - 1, j], Z), x[idx, j])$residuals
    cor(e, ex)
  })
  names(rho) <- cn
  structure(
    list(
      coefficients = beta, se = se, tstat = beta / se, vcov = Vb, Wald = W, df = k,
      p.value = 1 - stats::pchisq(W, k), gamma = cf[-seq_len(k)], rho_resid = rho,
      lags = q, robust = robust, n = length(idx), lm = fit
    ),
    class = "elliott_cf"
  )
}

# Eicker-White (HC0) covariance of an lm fit
hc0 <- function(fit) {
  X <- stats::model.matrix(fit)
  XtXi <- solve(crossprod(X))
  XtXi %*% crossprod(X * stats::residuals(fit)) %*% XtXi
}

#' @rdname elliott_cf
#' @param x an object of class "elliott_cf".
#' @param digits minimal number of significant digits.
#' @param ... unused.
#' @export
print.elliott_cf <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat("\nControl-function predictive regression (Elliott, 2011), ", x$lags, " covariate lag(s)\n\n", sep = "")
  tab <- cbind(Estimate = x$coefficients, "Std. Error" = x$se, "t value" = x$tstat,
               "Pr(>|t|)" = 2 * stats::pnorm(-abs(x$tstat)))
  printCoefmat(tab, digits = digits, P.values = TRUE, has.Pvalue = TRUE, signif.stars = FALSE)
  if (x$robust) cat("(Eicker-White standard errors)\n")
  cat("\nWald statistic: ", formatC(x$Wald, digits = digits), " on ", x$df, " DF, p-value ",
      format.pval(x$p.value, digits = digits), "\n", sep = "")
  cat("Remaining innovation correlation:", paste(names(x$rho_resid), formatC(x$rho_resid, digits = 3), collapse = ", "), "\n\n")
  invisible(x)
}

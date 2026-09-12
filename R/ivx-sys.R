#' Fitting Systems of IVX Predictive Regressions
#'
#' `ivx_sys` estimates a system of predictive regressions
#' \eqn{y_t = \mu + A x_{t-1} + \varepsilon_t} with an \eqn{m}-vector response
#' and \eqn{r} predictors of arbitrary persistence by IVX, and tests linear
#' restrictions on \eqn{A} with the IVX-Wald statistic of Kostakis,
#' Magdalinos and Stamatogiannis (2023), eqs (15) and (23), whose covariance
#' has the Kronecker form \eqn{(Z'X)^{-1} \otimes I_m} around
#' \eqn{Z(K)'Z(K) \otimes \hat\Sigma - n \bar z \bar z' \otimes \hat\Sigma_{FM}}.
#' Long horizons are handled as in [ivx()]. For \eqn{m = 1} the results
#' coincide with `ivx()`.
#'
#' @inheritParams ivx
#' @param formula a formula whose left-hand side is a matrix, e.g.
#' `cbind(y1, y2) ~ x1 + x2`.
#'
#' @return an object of class "ivx_sys" with the coefficient matrix
#' `coefficients` (responses in rows, predictors in columns), matching `se`,
#' `tstat` and `Wald_Ind` matrices, the joint Wald statistic `Wald_Joint`
#' (\eqn{\chi^2(mr)}), the per-equation Wald statistics `Wald_Eq`
#' (\eqn{\chi^2(r)}), and `vcov` of `vec(coefficients)` (column-major, names
#' `response:predictor`).
#'
#' @references Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023).
#' Taking stock of long-horizon predictability tests: Are factor returns
#' predictable? Journal of Econometrics, 237(2), 105380.
#'
#' @export
#' @examples
#' ivx_sys(cbind(Ret, DE) ~ DP + TBL, data = kms)
#'
#' summary(ivx_sys(cbind(Ret, DE) ~ DP + TBL, data = kms, horizon = 4))
ivx_sys <- function(formula, data, horizon, na.action, contrasts = NULL,
                    model = TRUE, x = FALSE, y = FALSE,
                    beta = 0.95, cz = 1, bandwidth = NULL, ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  if (missing(horizon)) horizon <- cl$horizon <- 1
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0) {
    warning("ivx estimation does not include an intercept by construction", call. = FALSE)
  }
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  if (!is.matrix(y)) y <- matrix(y, dimnames = list(NULL, deparse(formula[[2]])))
  x <- model.matrix(mt, mf, contrasts)
  z <- ivx_sys_fit(y, x, horizon = horizon, beta = beta, cz = cz, bandwidth = bandwidth, ...)
  class(z) <- "ivx_sys"
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model) z$model <- mf
  if (ret.x) z$x <- x
  if (ret.y) z$y <- y
  z
}

#' Fitter Function for Systems IVX Models
#'
#' Basic function called by `ivx_sys`. Should only be used directly by
#' experienced users.
#'
#' @param y numeric response matrix (observations in rows).
#' @inheritParams ivx_fit
#' @inheritParams ivx_sys
#' @export
#' @examples
#' ivx_sys_fit(as.matrix(kms[, c("Ret", "DE")]), as.matrix(kms[, c("DP", "TBL")]))
ivx_sys_fit <- function(y, x, horizon = 1, beta = 0.95, cz = 1, bandwidth = NULL, ...) {
  chkDots(...)
  y <- as.matrix(y)
  n <- NROW(x)
  if (is.null(n)) stop("'x' must be a matrix")
  if (NROW(y) != n) stop("incompatible dimensions")
  l <- NCOL(x)
  m <- NCOL(y)
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("x", 1L:l)
  rnames <- colnames(y)
  if (is.null(rnames)) rnames <- paste0("y", 1L:m)

  z <- ivx_sys_fit_cpp(y, x, horizon, beta, cz, bandwidth %||% -1L)

  dn <- list(rnames, cnames)
  coef <- z$A
  dimnames(coef) <- dimnames(z$se) <- dimnames(z$tstat) <- dn
  vnames <- as.vector(outer(rnames, cnames, paste, sep = ":"))
  dimnames(z$varcov) <- list(vnames, vnames)
  delta <- z$delta
  dimnames(delta) <- dn

  list(
    coefficients = coef,
    se = z$se,
    tstat = z$tstat,
    Wald_Ind = z$tstat^2,
    intercept = stats::setNames(drop(z$intercept), rnames),
    fitted = z$fitted,
    residuals = z$residuals,
    Wald_Joint = z$wald,
    Wald_Eq = stats::setNames(drop(z$wald_eq), rnames),
    rank = l,
    horizon = horizon,
    df.residuals = z$df.residuals,
    df = z$df,
    assign = attr(x, "assign"),
    cnames = cnames,
    rnames = rnames,
    AR = data.frame(Rn = drop(z$Rn), Rz = rep(z$Rz, l), row.names = cnames),
    delta = delta,
    vcov = z$varcov,
    tuning = list(beta = beta, cz = cz, bandwidth = z$bandwidth),
    ols = list(coefficients = z$ols$Aols, residuals = z$ols$residuals, Sigma = z$ols$Sigma)
  )
}

#' @rdname ivx_sys
#' @param x an object of class "ivx_sys".
#' @inheritParams stats::summary.lm
#' @export
print.ivx_sys <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients (responses in rows):\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname ivx_sys
#' @param object an object of class "ivx_sys".
#' @export
summary.ivx_sys <- function(object, ...) {
  z <- object
  est <- as.vector(z$coefficients)
  se <- as.vector(z$se)
  tt <- as.vector(z$tstat)
  coefs <- cbind(est, se, tt, tt^2, 1 - pchisq(tt^2, 1))
  dimnames(coefs) <- list(rownames(z$vcov), c("Estimate", "Std. Error", "t value", "Wald Ind", "Pr(> chi)"))
  ans <- list(
    call = z$call, coefficients = coefs, vcov = z$vcov, delta = z$delta,
    horizon = z$horizon, df = z$df, Wald_Joint = z$Wald_Joint,
    pv_waldjoint = 1 - pchisq(z$Wald_Joint, z$df),
    Wald_Eq = z$Wald_Eq, df_eq = length(z$cnames),
    pv_waldeq = 1 - pchisq(z$Wald_Eq, length(z$cnames))
  )
  class(ans) <- "summary.ivx_sys"
  ans
}

#' @export
print.summary.ivx_sys <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               cs.ind = 1:2, tst.ind = 3:4, signif.legend = TRUE, has.Pvalue = TRUE,
               P.values = TRUE, na.print = "NA", ...)
  cat("\nEquation Wald statistics on", x$df_eq, "DF:\n")
  print(data.frame(Wald = formatC(x$Wald_Eq, digits = digits),
                   "p-value" = format.pval(x$pv_waldeq, digits = digits),
                   row.names = names(x$Wald_Eq), check.names = FALSE))
  cat("\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
      "on", x$df, "DF, p-value", format.pval(x$pv_waldjoint, digits = digits), "\n\n")
  invisible(x)
}

#' @export
vcov.ivx_sys <- function(object, ...) object$vcov

#' @export
coef.ivx_sys <- function(object, ...) object$coefficients

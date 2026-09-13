#' Fitting IVX Quantile Predictive Regressions
#'
#' `ivx_qr` implements the IVX-QR predictability test of Lee (2016): the
#' \eqn{\tau}-quantile of the response is regressed on the IVX-filtered
#' predictors (Section 3.3, Proposition 3.2), giving a test of
#' \eqn{H_0: \beta_\tau = 0} with a standard chi-square limit whatever the
#' persistence of the predictors. Requires the \pkg{quantreg} package.
#'
#' The instrument is \eqn{z_t = (1 - c_z/n^\beta) z_{t-1} + \Delta x_t}. Lee
#' (2016) normalises \eqn{c_z = 5} and picks \eqn{\beta} from a look-up table
#' indexed by the estimated QR endogeneity
#' \eqn{\hat\rho(\tau) = -\mathrm{corr}(1\{\hat u_t < 0\}, \hat u_{x,t})}, which
#' is returned as `rho_tau` so the rule can be applied by the user; the default
#' `beta = 0.95` follows Kostakis et al. (2015). The sparsity
#' \eqn{f_u(0)} is estimated by a Gaussian kernel with Silverman's bandwidth
#' (footnote 4 of the paper).
#'
#' @inheritParams ivx
#' @param tau quantile level(s) in (0, 1). A vector fits one model per level.
#' @param ... further arguments passed to [quantreg::rq()].
#'
#' @return For a single `tau`, an object of class `c("ivx_qr", "ivx")` (so
#' `summary()`, `vcov()` etc. apply) with components `tau`, `rho_tau`,
#' `sparsity` (\eqn{\hat f_u(0)}) and `rq` (the underlying `quantreg::rq` fit).
#' For several `tau`, a list of such objects named by `tau`.
#'
#' @references Lee, J. H. (2016). Predictive quantile regression with persistent
#' covariates: IVX-QR approach. Journal of Econometrics, 192(1), 105-118.
#'
#' @export
#' @examples
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   summary(ivx_qr(Ret ~ DP, data = kms, tau = 0.5))
#'   ivx_qr(Ret ~ DP + TBL, data = kms, tau = c(0.1, 0.5, 0.9))
#' }
ivx_qr <- function(formula, data, tau = 0.5, beta = 0.95, cz = 5, na.action,
                   contrasts = NULL, model = TRUE, x = FALSE, y = FALSE, ...) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("package 'quantreg' is required for ivx_qr()", call. = FALSE)
  }
  if (any(tau <= 0 | tau >= 1)) stop("`tau` must be in (0, 1)", call. = FALSE)
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  fr <- ivx_frame(match.call(expand.dots = FALSE), parent.frame(), contrasts)

  fits <- lapply(tau, function(tt) {
    z <- ivx_qr_fit(fr$y, fr$x, tau = tt, beta = beta, cz = cz, ...)
    class(z) <- c("ivx_qr", "ivx")
    cl$tau <- tt
    ivx_finish(z, fr, cl, model, ret.x, ret.y)
  })
  if (length(tau) == 1) fits[[1]] else stats::setNames(fits, tau)
}

#' Fitter Function for IVX-QR Models
#'
#' Basic function called by `ivx_qr`. Should only be used directly by
#' experienced users.
#'
#' @inheritParams stats::lm.fit
#' @inheritParams ivx_qr
#' @export
#' @examples
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   ivx_qr_fit(kms$Ret, as.matrix(kms$DP), tau = 0.5)$Wald_Joint
#' }
ivx_qr_fit <- function(y, x, tau = 0.5, beta = 0.95, cz = 5, ...) {
  n <- NROW(x)
  l <- NCOL(x)
  if (is.null(n)) stop("'x' must be a matrix")
  if (NROW(y) != n) stop("incompatible dimensions")
  if (length(tau) != 1 || tau <= 0 || tau >= 1) stop("`tau` must be a single value in (0, 1)")
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("x", 1L:l)

  # IVX instrument z_t = rho z_{t-1} + dx_t (z_1 = 0); regress y_t on z_{t-1}, t = 2..n
  nn <- n - 1
  rho <- 1 - cz / nn^beta
  z <- apply(diff(x), 2, function(v) stats::filter(v, rho, "recursive"))
  zl <- rbind(0, matrix(z, ncol = l)[-nn, , drop = FALSE])  # z_{t-1}
  zl <- scale(zl, scale = FALSE)
  zl <- matrix(zl, ncol = l, dimnames = list(NULL, cnames))
  yt <- y[-1]

  fit <- quantreg::rq(yt ~ zl, tau = tau, ...)
  coef <- stats::setNames(coef(fit)[-1], cnames)
  u <- residuals(fit)

  # sparsity f_u(0): Gaussian kernel, Silverman's bandwidth
  h <- stats::bw.nrd0(u)
  f0 <- mean(stats::dnorm(u / h)) / h

  # Proposition 3.2: f0^2 / (tau (1 - tau)) * b' Z'Z b ~ chi^2(K)
  V <- tau * (1 - tau) / f0^2 * solve(crossprod(zl))
  dimnames(V) <- list(cnames, cnames)
  se <- sqrt(diag(V))
  tstat <- coef / se

  # QR endogeneity rho(tau) = -corr(1(u < 0), u_x), u_x from an AR(1) of each predictor
  ux <- sapply(seq_len(l), function(j) lm.fit(cbind(1, x[-n, j]), x[-1, j])$residuals)
  rho_tau <- stats::setNames(-drop(cor(as.numeric(u < 0), ux)), cnames)

  list(
    coefficients = coef,
    se = se,
    tstat = tstat,
    intercept = unname(coef(fit)[1]),
    fitted = drop(zl %*% coef),
    residuals = u,
    Wald_Joint = drop(crossprod(coef, solve(V, coef))),
    Wald_Ind = tstat^2,
    rank = l,
    horizon = 1,
    df.residuals = nn - l,
    df = l,
    assign = attr(x, "assign"),
    cnames = cnames,
    tau = tau,
    rho_tau = rho_tau,
    sparsity = f0,
    AR = data.frame(Rz = rep(rho, l), row.names = cnames),
    delta = matrix(rho_tau, 1, l, dimnames = list(NULL, cnames)),
    vcov = V,
    robust = FALSE,
    tuning = list(beta = beta, cz = cz, bandwidth = NA),
    data = list(y = yt, Z = zl),
    rq = fit
  )
}

#' Moving Block Bootstrap for IVX-QR
#'
#' Percentile confidence intervals and p-values for the IVX-QR coefficients
#' from the moving block bootstrap (MBB) of Fan and Lee (2019, Section 5).
#' Blocks of the pairs \eqn{(y_t, 	ilde z_{t-1})} are resampled and the
#' quantile regression of Lee (2016) is refitted on each bootstrap sample. This
#' avoids estimating the sparsity and the nuisance parameters that appear under
#' conditional heteroskedasticity, which is where the asymptotic IVX-QR test is
#' most distorted, especially in the tails.
#'
#' @param object an object of class "ivx_qr".
#' @param B number of bootstrap replications.
#' @param block block length; the default is \eqn{\lceil n^{1/4} \rceil} as in
#' the paper.
#' @param level confidence level of the percentile intervals.
#' @param seed optional integer seed.
#'
#' @return an object of class "ivx_qr_boot": a list with the estimates, the
#' percentile intervals `ci`, two-sided percentile p-values `p.value` for
#' \eqn{H_0: eta_j = 0}, and the bootstrap draws `boot`.
#'
#' @references Fan, R., & Lee, J. H. (2019). Predictive quantile regressions
#' under persistence and conditional heteroskedasticity. Journal of
#' Econometrics, 213(1), 261-280.
#'
#' @export
#' @examples
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   m <- ivx_qr(Ret ~ DP, data = kms, tau = 0.1)
#'   ivx_qr_boot(m, B = 199, seed = 1)
#' }
ivx_qr_boot <- function(object, B = 999, block = NULL, level = 0.95, seed = NULL) {
  if (!inherits(object, "ivx_qr")) stop("`object` must be an 'ivx_qr' fit", call. = FALSE)
  y <- object$data$y
  Z <- object$data$Z
  n <- length(y)
  b <- if (is.null(block)) ceiling(n^(1 / 4)) else as.integer(block)
  if (b < 1 || b > n) stop("`block` must be between 1 and the sample size", call. = FALSE)
  m <- ceiling(n / b)
  q <- n - b + 1
  if (!is.null(seed)) set.seed(seed)
  starts <- matrix(sample.int(q, m * B, replace = TRUE), B, m)
  offs <- seq_len(b) - 1L
  boot <- t(apply(starts, 1, function(st) {
    idx <- as.vector(outer(offs, st, "+"))            # m blocks of length b
    quantreg::rq.fit(cbind(1, Z[idx, , drop = FALSE]), y[idx], tau = object$tau)$coefficients[-1]
  }))
  boot <- matrix(boot, B, NCOL(Z), dimnames = list(NULL, object$cnames))
  a <- (1 - level) / 2
  ci <- t(apply(boot, 2, stats::quantile, probs = c(a, 1 - a)))
  # two-sided percentile p-value: smallest level at which 0 leaves the interval
  pv <- 2 * pmin(colMeans(boot <= 0), colMeans(boot >= 0))
  structure(
    list(
      call = object$call, tau = object$tau, B = B, block = b, level = level,
      coefficients = object$coefficients, ci = ci, p.value = pmin(pv, 1), boot = boot
    ),
    class = "ivx_qr_boot"
  )
}

#' @rdname ivx_qr_boot
#' @param x an object of class "ivx_qr_boot".
#' @param digits minimal number of significant digits.
#' @param ... unused.
#' @export
print.ivx_qr_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("
Call:
", paste(deparse(x$call), collapse = "
"), "

", sep = "")
  cat("IVX-QR at tau = ", x$tau, ", moving block bootstrap, B = ", x$B,
      ", block length ", x$block, "

", sep = "")
  tab <- cbind(Estimate = x$coefficients, x$ci, "Pr(|b| > 0)" = x$p.value)
  colnames(tab)[2:3] <- paste0(100 * c((1 - x$level) / 2, 1 - (1 - x$level) / 2), "%")
  cat("Coefficients (percentile intervals and p-values):
")
  printCoefmat(tab, digits = digits, cs.ind = 1:3, tst.ind = integer(),
               P.values = TRUE, has.Pvalue = TRUE, signif.stars = FALSE)
  cat("
")
  invisible(x)
}

#' @rdname ivx_qr
#' @param x an object of class "ivx_qr".
#' @inheritParams stats::summary.lm
#' @export
print.ivx_qr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("IVX-QR at tau = ", x$tau, "\n\nCoefficients:\n", sep = "")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname ivx_qr
#' @param object an object of class "ivx_qr".
#' @export
summary.ivx_qr <- function(object, ...) {
  ans <- summary.ivx(object, ...)
  ans$tau <- object$tau
  ans$rho_tau <- object$rho_tau
  class(ans) <- c("summary.ivx_qr", "summary.ivx")
  ans
}

#' @export
print.summary.ivx_qr <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("IVX-QR at tau = ", x$tau, "\n\nCoefficients:\n", sep = "")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               cs.ind = 1:2, tst.ind = 3:4, signif.legend = TRUE, has.Pvalue = TRUE,
               P.values = TRUE, na.print = "NA", ...)
  cat("\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
      "on", x$df, "DF, p-value", format.pval(x$pv_waldjoint, digits = digits))
  cat("\nQR endogeneity rho(tau):", paste(names(x$rho_tau), formatC(x$rho_tau, digits = 3), collapse = ", "))
  cat("\n\n")
  invisible(x)
}

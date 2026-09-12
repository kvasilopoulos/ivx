#' Fitting Residual-Augmented IVX Models
#'
#' `ivx_ra` implements the residual-augmented (bias-reduced) IVX estimator of
#' Demetrescu and Rodrigues (2022). An autoregression of order `p` is fitted to
#' the predictors, the predictive regression is augmented with its residuals
#' (in the spirit of Amihud and Hurvich, 2004), and the slope on the lagged
#' predictors is estimated by IVX. Inference uses the heteroskedasticity-robust
#' standard errors of the paper (eq. 9 and 14), which are valid whether the
#' predictors are stationary or near-integrated. Short horizon only.
#'
#' The autoregression of the predictors is fitted without an intercept and its
#' residuals are demeaned before augmentation, which is the paper's preferred
#' \eqn{\tilde t_{ivx}^{\mu_0}} statistic (Sections 4-5); the standard errors
#' include the finite-sample correction of Kostakis et al. (2015) as in the
#' paper's simulations.
#'
#' @inheritParams ivx
#' @param ar order of the autoregression fitted to the predictors: `"auto"`
#' selects it by `ar_ic` in levels (as recommended in the paper), or a positive
#' integer for a fixed order.
#' @param ar_ic information criterion for `ar = "auto"`.
#' @param ar_max maximum order considered when `ar = "auto"`.
#'
#' @return an object of class `c("ivx_ra", "ivx")`; the usual `ivx` methods
#' (`summary`, `vcov`, `delta`, ...) apply. Additional components: `gamma`
#' (coefficients on the augmentation residuals) and `ar_order`.
#'
#' @references Demetrescu, M., & Rodrigues, P. M. M. (2022). Residual-augmented
#' IVX predictive regression. Journal of Econometrics, 227(2), 429-460.
#' @references Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
#' reduced-bias estimation method. Journal of Financial and Quantitative
#' Analysis, 39(4), 813-841.
#'
#' @export
#' @examples
#' ivx_ra(Ret ~ DP, data = kms)
#'
#' summary(ivx_ra(Ret ~ DP + TBL, data = kms, ar = 2))
ivx_ra <- function(formula, data, ar = "auto", ar_ic = c("aic", "bic"), ar_max = 5,
                   beta = 0.95, cz = 1, na.action, contrasts = NULL,
                   model = TRUE, x = FALSE, y = FALSE, ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  ar_ic <- match.arg(ar_ic)
  if (!identical(ar, "auto") && !(is.numeric(ar) && length(ar) == 1 && ar >= 1 && ar == trunc(ar))) {
    stop("`ar` should be either 'auto' or a positive integer.", call. = FALSE)
  }
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
  if (is.matrix(y)) {
    stop("multivariate model is not available", call. = FALSE)
  }
  x <- model.matrix(mt, mf, contrasts)
  z <- ivx_ra_fit(y, x, ar = ar, ar_ic = ar_ic, ar_max = ar_max, beta = beta, cz = cz, ...)
  class(z) <- c("ivx_ra", "ivx")
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

#' Fitter Function for Residual-Augmented IVX Models
#'
#' Basic function called by `ivx_ra` to fit the model. Should only be used
#' directly by experienced users.
#'
#' @inheritParams stats::lm.fit
#' @inheritParams ivx_ra
#' @param ... currently unused.
#' @export
#' @examples
#' ivx_ra_fit(kms$Ret, as.matrix(kms$DP))
ivx_ra_fit <- function(y, x, ar = "auto", ar_ic = "aic", ar_max = 5, beta = 0.95, cz = 1, ...) {
  chkDots(...)
  n <- NROW(x)
  l <- NCOL(x)
  if (is.null(n)) stop("'x' must be a matrix")
  if (NROW(y) != n) stop("incompatible dimensions")
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("x", 1L:l)

  # Step 1: AR(p) / VAR(p) in levels without intercept (mu_0); residuals eps_t, t = p+1..n
  va <- if (identical(ar, "auto")) {
    var_ols(x, ar_max, ic = ar_ic, intercept = FALSE)
  } else {
    var_ols(x, ar, fixed = TRUE, intercept = FALSE)
  }
  p <- va$q
  # demeaned, i.e. regression (6) with an intercept; without this the term
  # gamma * mean(eps) is not removed and the bias reduction fails
  eps <- scale(va$resid, scale = FALSE)
  eps <- matrix(eps, ncol = l)
  idx <- (p + 1):n
  nn <- length(idx)

  # IVX instrument as in Kostakis et al. (2015): z_1 = 0, z_t = rho z_{t-1} + dx_t
  rho <- 1 - cz / (n - 1)^beta
  zf <- rbind(0, apply(diff(x), 2, function(v) stats::filter(v, rho, "recursive")))
  zl <- matrix(zf[idx - 1, ], ncol = l)                 # z_{t-1}, not demeaned
  yd <- y[idx] - mean(y[idx])
  xd <- scale(x[idx - 1, , drop = FALSE], scale = FALSE)  # demeaned x_{t-1}
  xd <- matrix(xd, ncol = l)

  # Step 2: OLS of demeaned y on the AR residuals, partial them out
  gamma <- lm.fit(eps, yd)$coefficients
  ytil <- drop(yd - eps %*% gamma)

  # Step 3: IVX of ytil on x_{t-1}
  B <- crossprod(zl, xd)
  coef <- drop(solve(B, crossprod(zl, ytil)))

  # robust covariance, eq. (9) / (14): residuals from the OLS augmented regression
  ols <- lm.fit(cbind(xd, eps), yd)
  e <- ols$residuals
  Xp <- scale(embed(x, p + 1)[, -seq_len(l), drop = FALSE], scale = FALSE)  # x_{t-1..t-p}, demeaned
  Xp <- matrix(Xp, ncol = p * l)
  ge <- drop(eps %*% gamma)                                # gamma' eps_t
  Hxx <- crossprod(Xp)
  Hzx <- crossprod(zl, Xp)
  corr <- Hzx %*% solve(Hxx, crossprod(Xp * ge)) %*% solve(Hxx, t(Hzx))
  # plus the Kostakis et al. (2015) intercept correction, as in the paper's Section 4
  M <- crossprod(zl * e) + corr - nn * tcrossprod(colMeans(zl)) * kms_fm(e, eps)
  Binv <- solve(B)
  V <- Binv %*% M %*% t(Binv)

  se <- sqrt(diag(V))
  tstat <- coef / se
  names(coef) <- names(se) <- names(tstat) <- cnames
  dimnames(V) <- list(cnames, cnames)
  fitted <- drop(xd %*% coef)

  u_ols <- lm.fit(xd, yd)$residuals
  delta <- matrix(cor(u_ols, eps), 1, l, dimnames = list(NULL, cnames))

  list(
    coefficients = coef,
    se = se,
    tstat = tstat,
    gamma = stats::setNames(drop(gamma), cnames),
    intercept = mean(y[idx]) - sum(colMeans(x[idx - 1, , drop = FALSE]) * coef),
    fitted = fitted,
    residuals = ytil - fitted,
    Wald_Joint = drop(crossprod(coef, solve(V, coef))),
    Wald_Ind = tstat^2,
    rank = l,
    horizon = 1,
    df.residuals = nn - l,
    df = l,
    assign = attr(x, "assign"),
    cnames = cnames,
    ar_order = p,
    ar_method = if (identical(ar, "auto")) ar_ic else "fixed",
    AR = data.frame(Rz = rep(rho, l), row.names = cnames),
    delta = delta,
    vcov = V,
    robust = TRUE,
    tuning = list(beta = beta, cz = cz, bandwidth = NA)
  )
}

#' @rdname ivx_ra
#' @param x an object of class "ivx_ra".
#' @inheritParams stats::summary.lm
#' @export
print.ivx_ra <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Residual-augmented IVX, AR order p = ", x$ar_order,
      if (x$ar_method == "fixed") " (fixed)" else paste0(" (", x$ar_method, ")"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

# sigma_e^2 - Omega_eu' Omega_uu^-1 Omega_eu with Bartlett/Newey-West long-run
# estimates (bandwidth n^(1/3)), as in ivx_fit_cpp; e: n-vector, u: n x l matrix
kms_fm <- function(e, u) {
  n <- length(e)
  m <- floor(n^(1 / 3))
  w <- 1 - seq_len(m) / (1 + m)
  cov_uu <- crossprod(u) / n
  cov_eu <- crossprod(u, e) / n
  uu <- 0
  eu <- 0
  for (h in seq_len(m)) {
    uu <- uu + w[h] * crossprod(u[-seq_len(h), , drop = FALSE], u[seq_len(n - h), , drop = FALSE]) / n
    eu <- eu + w[h] * crossprod(u[-seq_len(h), , drop = FALSE], e[seq_len(n - h)]) / n
  }
  omega_uu <- cov_uu + uu + t(uu)
  omega_eu <- cov_eu + eu
  drop(mean(e^2) - crossprod(omega_eu, solve(omega_uu, omega_eu)))
}

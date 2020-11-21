# These functions have been refactored from forecast::auto.arima

auto_ar <- function(x, max.p = 5, ic = "aic", ...) {
  p <- max.p
  bestfit <- arima2(x, order = c(p, 0, 0), include.mean = FALSE)
  while (p >= 0) {
    fit <- arima2(x, order = c(p, 0, 0), include.mean = FALSE)
    if (fit[[ic]] < bestfit[[ic]]) {
      bestfit <- fit
    }
    p <- p - 1
  }
  bestfit
}

#' @importFrom stats arima frequency
arima2 <- function(x, order = c(0, 0, 0), ic = "aic", ...) {
  fit <- arima(x = x, order = order, ...)
  npar <- length(fit$coef[fit$mask]) + 1
  n <- NROW(x)
  m <- frequency(x)
  nstar <- n - m
  # dots <- list(...)
  # if (is.null(dots$xreg)) {
  #   nxreg <- 0
  # } else {
  #   nxreg <- ncol(as.matrix(xreg))
  # }
  if (!is.na(fit$aic)) {
    fit$bic <- fit$aic + npar * (log(nstar) - 2)
    fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
    fit$ic <- switch(ic, bic = fit$bic, aic = fit$aic, aicc = fit$aicc)
  }
  else {
    fit$aic <- fit$bic <- fit$aicc <- fit$ic <- Inf
  }
  fit
}

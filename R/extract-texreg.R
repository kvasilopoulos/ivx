


#' \code{\link{extract}} method for \code{ivx} objects
#'
#' \code{\link{extract}} method for \code{Arima} objects created by the
#' \code{\link[stats]{arima}} function in the \pkg{stats} package.
#'
#' @param model A statistical model object.
#' @param include.pvalues Report p-values?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.loglik Report the log likelihood in the GOF block?
#' @param include.nobs Report the number of observations in the GOF block?
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#'
#' @export
extract.ivx <- function(model, include.wald = TRUE,
                        inlude.nobs = TRUE,
                        include.rsquared = FALSE,
                        include.adjrs = FALSE) {
  s <- summary(model)

  names <- rownames(s$coef)
  co <- s$coef[, 1]
  wald <- s$coef[, 2]
  pval <- s$coef[, 3]

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.wald == TRUE) {
    wl <- s$Wald_Joint
    gof <- c(gof, wl)
    gof.names <- c(gof.names, "Joint Wald$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared == TRUE) {
    rs <- s$r.squared
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    adj <- s$adj.r.squared
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if(inlude.nobs == TRUE) {
    n <- NROW(model$residuals)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- texreg::createTexreg(
    coef.names = names,
    coef = co,
    se = wald,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' @export
extract.ivx_ar <- extract.ivx

.onLoad <- function(libname, pkgname) {
  if (suppressWarnings(requireNamespace("texreg", quietly = TRUE))) {
    setGeneric("extract", function(model, ...) standardGeneric("extract"),
               package = "texreg"
    )
    setMethod("extract", signature = className("ivx", "ivx"),
              definition = extract.ivx)
  }
}


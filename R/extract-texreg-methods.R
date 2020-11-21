#' \code{extract} method for \code{ivx} objects
#'
#'
#' @param model A statistical model object.
#' @param include.wald Report the Wald statistic.
#' @param include.nobs Report the number of observations in the GOF block?
#' @param include.aic Report Akaike's Information Criterion (AIC) in the GOF
#'   block?
#' @param include.bic Report the Bayesian Information Criterion (BIC) in the GOF
#'   block?
#' @param include.rsquared Report the R-squared.
#' @param include.adjrs Report the Adjusted R-squared.
#' @param ... Custom parameters, which are handed over to subroutines. Currently
#'   not in use.
#'
#' @importFrom stats AIC BIC
#' @export
extract.ivx <- function(model,
                        include.wald = TRUE,
                        include.nobs = TRUE,
                        include.aic = FALSE,
                        include.bic = FALSE,
                        include.rsquared = FALSE,
                        include.adjrs = FALSE,
                        ...) {
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
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, BIC(model))
    gof.names <- c(gof.names, "BIC")
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
  if(include.nobs == TRUE) {
    n <- NROW(model$residuals)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  texreg::createTexreg(
    coef.names = names,
    coef = co,
    se = wald,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
}


#' @rdname extract.ivx
#' @importFrom methods setGeneric setMethod className
#' @export
extract.ivx_ar <- extract.ivx

.onLoad <- function(libname, pkgname) {
  if (suppressWarnings(requireNamespace("texreg", quietly = TRUE))) {
    methods::setGeneric("extract", function(model, ...) standardGeneric("extract"),
               package = "texreg"
    )
    methods::setMethod("extract", signature = methods::className("ivx", "ivx"),
              definition = extract.ivx)
  }
}


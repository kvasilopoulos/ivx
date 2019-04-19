#' Fitting IVX Models
#'
#' @inheritParams stats::lm
#' @param horizon is the horizon (default horizon=1 corresponds to a
#' short-horizon regression)
#'
#' @return ivx returns an object of class "ivx".
#'
#' An object of class "ivx" is a list containing the following components:
#' \item{coefficients_ols}{coefficients}
#' \item{coefficients_ivx}{coefficients}
#' \item{call}{call}
#' \item{terms}{terms}
#' \item{model}{ model}
#'
#' @references Phillips, P. C., & Magdalinos, T. (2009). Econometric inference
#' in the vicinity of unity. Singapore Management University, CoFie Working
#' Paper, 7.
#' @references Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2014).
#' Robust econometric inference for stock return predictability. The Review of
#' Financial Studies, 28(5), 1506-1553.
#' @references Phillips, P. C., & Lee, J. H. (2013). Predictive regression under
#' various degrees of persistence and robust long-horizon regression. Journal
#' of Econometrics, 177(2), 250-264.
#'
#' @export
#'
#' @importFrom stats .getXlevels coef coefficients cor lm model.matrix pf
#' model.offset model.response pchisq qnorm residuals symnum
ivx <- function(formula, data, horizon, subset, na.action,
                 model = TRUE, contrasts = NULL, offset, ...)
{

  cl <- match.call()

  if (missing(horizon)) horizon <- cl$horizon <- 1

  ## keep only the arguments which should go into the model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "horizon","subset", "na.action",
               "offset"), names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame) # was as.name("model.frame"), but
  ##    need "stats:: ..." for non-standard evaluation
  mf["horizon"] <- NULL
  # mf$formula
  mf <- eval.parent(mf)
  # if (method == "model.frame") return(mf)

  ## 1) allow model.frame to update the terms object before saving it.
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0) {
    warning("ivx estimation does not include an intercept by construction",
            call. = FALSE)
  }
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")

  ## 2) retrieve the weights and offset from the model frame so
  ## they can be functions of columns in arg data.
  # w <- model.weights(mf)
  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)
  ## if any subsetting is done, retrieve the "contrasts" attribute here.

  z <- ivx_fit(y, x, h = horizon, offset, ...) # offset = offset,
  class(z) <- "ivx"

  ## 3) return the na.action info
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset

  ## 4) return the contrasts used in fitting: possibly as saved earlier.
  z$contrasts <- attr(x, "contrasts")

  ## 5) return the levelsets for factors in the formula
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)  z$model <- mf
  z
}



#' @rdname ivx
#' @export
ivx_fit <- function(y, x, h = 1, offset = NULL, ...) {

  n <- NROW(x)
  p <- NCOL(x)

  if (is.null(n)) stop("'x' must be a matrix")
  if (n == 0L)  stop("0 (non-NA) cases")
  if (p == 0L) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, rank = 0, df.residual = length(y)))
  }

  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1)
    y <- drop(y)
  if (!is.null(offset))
    y <- y - offset
  if (NROW(y) != n)  stop("incompatible dimensions")

  chkDots(...)

  z <- ivx_fit_cpp(y, x, h)

  cnames <- colnames(x)
  coef <- drop(z$Aivx)
  coef_ols <- drop(z$Aols)

  if (is.null(cnames)) cnames <- paste0("x", 1L:p)
  # nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
  z$coefficients <- coef
  z$coefficients_ols <- coef_ols
  # r1 <- y - z$residuals

  if (!is.null(offset)) r1 <- r1 + offset

  if (is.matrix(y)) {
    dimnames(coef) <- list(cnames, colnames(y))
    dimnames(coef_ols) <- list(c("Intercept", cnames), colnames(y))
    # dimnames(z$effects) <- list(nmeffects, colnames(y))
  }else{
    names(coef) <- cnames
    names(coef_ols) <- c("Intercept", cnames)
  }

  output <-
    structure(
      list(coefficients =  coef,
           fitted = drop(z$fitted),
           residuals = drop(z$residuals),
           Wald_Joint = z$wivx,
           Wald_Ind = z$wivxind,
           horizon = h,
           df = z$df,
           cnames = cnames,
           AR = data.frame(Rn = z$Rn,
                           Rz = z$Rz,
                           row.names = cnames),
           delta = z$delta,
           vcov = z$varcov,
           coefficients_ols = coef_ols,
           tstat_ols = z$tstat_ols
      )
    )
  output
}

#' @rdname ivx
#' @export
print.ivx <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  res <- x$coefficients
  if (length(res)) {
    cat("Coefficients:\n")
    print.default(format(res, digits = digits), print.gap = 2L, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  invisible(x)
}

#' Summarizing IVX Model Fits
#'
#' summary method for class "ivx".
#'
#' @param object  object of class "ivx", usually, a result of a call to ivx.
#' @inheritParams stats::summary.lm
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @importFrom stats printCoefmat
#' @importFrom stats pt
summary.ivx <- function(object,  ...) {
  z <- object

  if (is.null(z$terms))
    stop("invalid 'ivx' object: no 'terms' components")
  if (!inherits(object, "ivx"))
    stop("calling summary.ivx(<fake-ivx-object>) ...")

  ans <- z[c("call", "terms")]

  ans$aliased <- is.na(z$coefficients)

  p_value_ivx <- 1 - pchisq(z$Wald_Ind, 1)
  ans$coefficients <- cbind(z$coefficients, z$Wald_Ind, p_value_ivx)
  dimnames(ans$coefficients) <- list(z$cnames, c("Estimate", "Wald Ind", "Pr(> chi)"))


  ans$vcov <- z$vcov
  dimnames(ans$vcov) <- dimnames(ans$coefficients)[c(1, 1)]

  ans$delta <- z$delta
  colnames(ans$delta) <- z$cnames

  ans$residuals <- z$residuals
  ans$fitted <- z$fitted

  ans$horizon <- z$horizon
  ans$Wald_Joint <- z$Wald_Joint
  ans$pv_waldjoint <- 1 - pchisq(z$Wald_Joint, z$df[1]);
  ans$df <- z$df

  if (is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.ivx"

  ans
}

#' @inheritParams stats::summary.lm
#' @rdname summary.ivx
#' @export
print.summary.ivx <- function(x,
                              digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars"),
                              ...){

  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (length(x$aliased) == 0L) {
    cat("No coefficients\n")
  }else{

    coefs_ols <- x$coefficients_ols
    coefs <- x$coefficients
    aliased <- x$aliased

    if (!is.null(aliased) && any(aliased)) {
      cn <- names(aliased)
      civx <- x$coefficients_ols
      coefs <- matrix(NA, NROW(civx), 5, dimnames = list(cn , colnames(civx)))
      coefs[!aliased, ] <- civx
    }

    cat("Coefficients:\n")

    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 signif.legend = TRUE, has.Pvalue = TRUE, P.values = TRUE,
                 na.print = "NA", ...)

    cat("\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
        "on", x$df[1], "DF, p-value",
        format.pval(x$pv_waldjoint, digits = digits))
  }

 invisible(x)
}

#' Calculate the delta coefficient
#'
#' Calculate the correlation coefficient between the residuals of regression between
#' the predictive and the the autoregressive regression.
#'
#' @param object class "ivx"
#'
#' @return This should have row and column names corresponding to the parameter
#' names given by the coef method.
#'
#' @export
delta <- function(object) {

  if (!inherits(object, c("ivx", "summary.ivx"))) {
    stop("Wrong object", call. = FALSE)
  }
  drop(object[["delta"]])

}

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' @inheritParams stats::vcov.lm
#' @param complete logical indicating if the full variance-covariance matrix
#' should be returned also in case of an over-determined system where some
#' coefficients are undefined and coef(.) contains NAs correspondingly. When
#' complete = TRUE, vcov() is compatible with coef() also in this singular case.
#' @param ... additional arguments for method functions.
#'
#' @return A matrix of the estimated covariances between the parameter estimates
#' in the linear or non-linear predictor of the model. This should have row and
#' column names corresponding to the parameter names given by the coef method.
#'
#' @export
vcov.ivx <- function(object, complete = TRUE, ...) {
  vcov.summary.ivx(summary.ivx(object), complete = complete, ...)
}

#' @rdname vcov.ivx
#' @export
vcov.summary.ivx <- function(object, complete = TRUE, ...) {
  stats::.vcov.aliased(object$aliased, object$vcov, complete = complete)
}

#' @export
residuals.ivx <- function(object, ...) {
  residuals.summary.ivx(summary.ivx(object), ...)
}
#' @importFrom stats naresid
#' @export
residuals.summary.ivx <- function(object, ...) {
  r <- object$residuals
  stats::naresid(object$na.action, r)
}



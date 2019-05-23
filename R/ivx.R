#' Fitting IVX Models
#'
#' ivx fits predictive regression models. The method allows standard
#' chi-square testing for regressors with different degrees of persistence,
#' from stationary to mildly explosive, and can be used for both short-
#' and long-horizon predictive regressions.
#'
#' @inheritParams stats::lm
#' @param horizon is the horizon (default horizon = 1 corresponds to a
#' short-horizon regression)
#'
#' @return an object of class "ivx".
#'
#' @references Magdalinos, T., & Phillips, P. (2009). Limit Theory for Cointegrated
#' Systems with Moderately Integrated and Moderately Explosive Regressors.
#' Econometric Theory, 25(2), 482-526.
#' @references Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2014).
#' Robust econometric inference for stock return predictability. The Review of
#' Financial Studies, 28(5), 1506-1553.
#'
#' @aliases ivx
#'
#' @importFrom stats .getXlevels coef coefficients cor lm model.matrix pf
#' model.offset model.response pchisq qnorm residuals symnum
#'
#' @export
#' @examples
#'
#' # Univariate
#' ivx(Ret ~ LTY, data = monthly)
#'
#' # Multivariate
#' ivx(Ret ~ LTY + TBL , data = monthly)
#'
#' # Longer horizon
#' ivx(Ret ~ LTY + TBL, data = monthly, horizon = 4)
ivx <- function(formula, data, horizon, na.action,
                contrasts = NULL, offset, ...)
{

  cl <- match.call()

  if (missing(horizon)) horizon <- cl$horizon <- 1

  ## keep only the arguments which should go into the model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "horizon", "na.action",
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

  z <- ivx_fit(y, x, horizon = horizon, offset, ...) # offset = offset,
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
  z
}



#' Fitter Functions for IVX Models
#'
#' Basic function called by `ivx` to fit predictive models.
#' These should only be used directly by experienced users.
#'
#' @inheritParams stats::lm.fit
#' @inheritParams ivx
#' @export
#' @examples
#' ivx_fit(monthly$Ret, as.matrix(monthly$LTY))
ivx_fit <- function(y, x, horizon = 1, offset = NULL, ...) {

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

  z <- ivx_fit_cpp(y, x, horizon)

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
           horizon = horizon,
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
#' @inheritParams stats::summary.lm
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
#' @examples
#' mod <- ivx(Ret ~ LTY, data = monthly)
#'
#' summary(mod)
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
  ans$pv_waldjoint <- 1 - pchisq(z$Wald_Joint, z$df[1])
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
#' Computes the long-run correlation coefficient between the residuals of the
#' predictive regression and the autoregressive model for the regressor.
#'
#' @param object on object of class "ivx"
#'
#' @return A vector of the estimated correlation coefficients. This should have
#' row and column names corresponding to the parameter names given by the coef method.
#'
#' @export
#' @examples
#' mod <- ivx(Ret ~ LTY, data = monthly)
#'
#' delta(mod)
delta <- function(object) {

  if (!inherits(object, c("ivx", "summary.ivx"))) {
    stop("Wrong object", call. = FALSE)
  }
  drop(object[["delta"]])

}

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' @param object a fitted ivx and summary.ivx object.
#' @param complete logical indicating if the full variance-covariance matrix
#' should be returned. When complete = TRUE, vcov() is compatible with coef().
#' @param ... additional arguments for method functions.
#'
#' @return A matrix of the estimated covariances between the parameter estimates
#' of the model. This should have row and column names corresponding to the
#' parameter names given by the coef method.
#'
#' @export
#' @examples
#' mod <- ivx(Ret ~ LTY, data = monthly)
#'
#' vcov(mod)
vcov.ivx <- function(object, complete = TRUE, ...) {
  vcov.summary.ivx(summary.ivx(object), complete = complete, ...)
}

#' @rdname vcov.ivx
#' @export
vcov.summary.ivx <- function(object, complete = TRUE, ...) {
  stats::.vcov.aliased(object$aliased, object$vcov, complete = complete)
}



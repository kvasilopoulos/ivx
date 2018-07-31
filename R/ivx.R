#' IVX estimation
#'
#' Date-filtering procedure that removes endogeneity.
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
#' @importFrom pracma pinv
#' @importFrom stats .getXlevels coef coefficients cor lm model.matrix
#' model.offset model.response pchisq qnorm residuals symnum
ivx <- function(formula, data, horizon, subset, na.action,
                 model = TRUE, contrasts = NULL, offset, ...)
{

  cl <- match.call()
  # attr(formula, "intercept") <- 0

  if (missing(horizon)) {
    horizon = 1
    cl$horizon <- 1
  }

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
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")

  ## 2) retrieve the weights and offset from the model frame so
  ## they can be functions of columns in arg data.
  # w <- model.weights(mf)
  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)
  ## if any subsetting is done, retrieve the "contrasts" attribute here.

  z <- ivx.fit(y, x, h = horizon, offset, ...) # offset = offset,
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

#' @export
print.ivx <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  res <- x$coefficients_ivx

  if (length(res)) {
    cat("Coefficients:\n")
    print.default(format(res, digits = digits),
                  print.gap = 2L, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  invisible(x)
  cat("\nNote: ivx estimation does not include an intercept.")
}

#' @export
#' @importFrom stats printCoefmat
summary.ivx <- function(object, ...) {
  z <- object

  if (is.null(z$terms))
    stop("invalid 'ivx' object: no 'terms' components")
  if (!inherits(object, "ivx"))
    stop("calling summary.ivx(<fake-ivx-object>) ...")

  ans <- z[c("call", "terms")]
  ans$coefficients_ivx <- cbind(z$coefficients_ivx, z$WivxInd, z$WivxInd_pvalue,
                                z$delta, z$ar)
  dimnames(ans$coefficients_ivx) <- list(z$cnames, c("Estimate", "Wald Ind",
                                       "Pr(> chi)", "delta", "Rn"))#\u03c7
  ans$coefficients_ols <- z$coefficients_ols
  ans$aliased <- is.na(z$coefficients_ivx)
  ans$horizon <- z$horizon
  ans$Wivx <- z$Wivx
  ans$Wivx_pvalue <- z$Wivx_pvalue
  ans$sigma <- z$varcovIVX
  if (is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.ivx"
  #print(ans$coefficients_ivx, digits = max(3L, getOption("digits") - 3L))
  ans
}

#' @export
print.summary.ivx <- function(x, digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars"),
                              ...){
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(x$aliased) == 0L) {
    cat("No coefficients\n")
  }else{
    cat("Coefficients:\n")
    coefs <- x$coefficients_ivx
    aliased <- x$aliased

    if (!is.null(aliased) && any(aliased)) {
      cn <- names(aliased)
      civx <- x$coefficients_ivx
      coefs <- matrix(NA, NROW(civx), 5,dimnames = list(cn , colnames(civx)))
      coefs[!aliased, ] <- civx
    }

    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  }

 cat("\nJoint Wald statistic", format(x$Wivx, digits = digits),
     ", p-value", format(x$Wivx_pvalue, digits = digits))

 # cat("\n\n Coefficients OLS:\n")
 # printCoefmat(x$coefficients_ols, digits = digits, signif.stars = signif.stars,
 #              na.print = "NA", ...)
 invisible(x)
}


# coef.ivx <- function(x, ols = FALSE) {
#   if (ols) {
#     x$coefficients_ols
#   }else{
#     x$coefficients_ivx
#   }
# }

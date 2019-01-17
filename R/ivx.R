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
  res <- x$coefficients
  if (length(res)) {
    cat("Coefficients:\n")
    print.default(format(res, digits = digits), print.gap = 2L, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  invisible(x)
}

#' @export
#' @importFrom stats printCoefmat
summary.ivx <- function(object, ols = FALSE, ...) {
  z <- object

  if (is.null(z$terms))
    stop("invalid 'ivx' object: no 'terms' components")
  if (!inherits(object, "ivx"))
    stop("calling summary.ivx(<fake-ivx-object>) ...")

  ans <- z[c("call", "terms")]
  ans$coefficients <- cbind(z$coefficients, z$Wald_Ind) #, z$WivxInd_pvalue
  ans$delta <- z$delta
  # z$delta, z$ar
  dimnames(ans$coefficients) <- list(z$cnames, c("Estimate", "Wald Ind")) #,"Pr(> chi)")
                                         #\u03c7 #, "delta", "Rn")
  ans$coefficients_ols <- z$coefficients_ols
  # ans$fstatistic <- z$fstatistic
  ans$aliased <- is.na(z$coefficients)
  ans$horizon <- z$horizon
  ans$Wivx <- z$Wald_Joint

  ans$df <- z$df
  ans$pv_ivx <- 1 - pchisq(z$Wald_Ind, z$df);
  ans$pv_wivxind <- 1 - pchisq(z$Wald_Joint^2, 1);

  # ans$sigma <- z$varcovIVX

  if (is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.ivx"
  attr(ans, "ols") <- ols
  ans
}

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

    coefs <- x$coefficients_ols
    aliased <- x$aliased

    if (!is.null(aliased) && any(aliased)) {
      cn <- names(aliased)
      civx <- x$coefficients_ols
      coefs <- matrix(NA, NROW(civx), 5, dimnames = list(cn , colnames(civx)))
      coefs[!aliased, ] <- civx
    }

    if (attr(x, "ols")) {
      cat("\n\nCoefficients OLS:\n")

      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                   signif.legend = TRUE, na.print = "NA", ...)
      if (!is.null(x$fstatistic)) {
        cat("\nF-statistic:", formatC(x$fstatistic[1L], digits = digits),
            "on", x$fstatistic[2L], "and",
            x$fstatistic[3L], "DF,  p-value:",
            format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                           x$fstatistic[3L], lower.tail = FALSE),
                        digits = digits))
      }
    }else{

      coefs <- x$coefficients
      aliased <- x$aliased

      if (!is.null(aliased) && any(aliased)) {
        cn <- names(aliased)
        civx <- x$coefficients
        coefs <- matrix(NA, NROW(civx), 5, dimnames = list(cn , colnames(civx)))
        coefs[!aliased, ] <- civx
      }

      cat("Coefficients:\n")
      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                   if (attr(x, "ols")) signif.legend = TRUE else signif.legend = FALSE,
                   has.Pvalue = TRUE, P.values = TRUE, na.print = "NA", ...)

      cat("\nJoint Wald statistic: ", formatC(x$wivx, digits = digits),
          "on", x$df, "DF, p-value",
          format.pval(x$pv_wivx, digits = digits))
    }

  }

 invisible(x)
}

#' Finding the delta coefficient
#'
#' @param object class "ivx"
#' @param mat matrix format
#' @param diag include diagonal elements
#'
#' @importFrom gdata upperTriangle
#' @export
delta <- function(object, mat = F, diag = FALSE) {

  if (!inherits(object, c("ivx", "summary.ivx"))) {
    stop("Wrong object", call. = F)
  }
  delta <- object$delta
  if (mat) {
    upperTriangle(delta, diag = diag) <- NA
  }else{
    delta <- drop(delta[-1, 1])
  }
  delta
}

varcov <- function(object) {
  if (!inherits(x, c("ivx", "summary.ivx"))) stop("Wrong object", call. = F)
  varcov <- object$varcov

}

coef.ivx <- function(x, ols = FALSE) {
  if (ols) {
    x$coefficients
  }else{
    x$coefficients_ivx
  }
}

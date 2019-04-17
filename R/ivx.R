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


ivx_fit <- function(y, x, h = 1, offset = NULL, ...) {

  if (is.null(n <- NROW(x))) stop("'x' must be a matrix")
  if (n == 0L)  stop("0 (non-NA) cases")
  p <- NCOL(x)
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
           residuals = z$residuals,
           Wald_Joint = z$wivx,
           Wald_Ind = z$wivxind,
           horizon = h,
           df = z$df,
           cnames = cnames,
           AR = data.frame(Rn = z$Rn,
                           Rz = z$Rz,
                           row.names = cnames),
           delta = z$delta,
           varcov = z$varcov,
           coefficients_ols = coef_ols,
           tstat_ols = z$tstat_ols
      )
    )
  output
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
#' @importFrom stats pt
summary.ivx <- function(object, ols = FALSE, ...) {
  z <- object

  if (is.null(z$terms))
    stop("invalid 'ivx' object: no 'terms' components")
  if (!inherits(object, "ivx"))
    stop("calling summary.ivx(<fake-ivx-object>) ...")

  ans <- z[c("call", "terms")]

  ans$aliased <- if (ols) is.na(z$coefficients_ols) else is.na(z$coefficients)

  ans$delta <- z$delta

  # OLS
  p_value_ols <- 2*(1 - pt(abs(z$tstat_ols), z$df[2]))
  ans$coefficients_ols <- cbind(z$coefficients_ols, z$tstat_ols, p_value_ols)
  colnames(ans$coefficients_ols) <- c("Estimate", "t value", "Pr(> |t|)")

  # IVX
  p_value_ivx <- 1 - pchisq(z$Wald_Ind, 1)
  ans$coefficients <- cbind(z$coefficients, z$Wald_Ind, p_value_ivx)
  dimnames(ans$coefficients) <- list(z$cnames, c("Estimate", "Wald Ind", "Pr(> chi)"))


  ans$horizon <- z$horizon
  ans$Wald_Joint <- z$Wald_Joint
  ans$pv_waldjoint <- 1 - pchisq(z$Wald_Joint, z$df[1]);
  ans$sigma <- z$varcov
  ans$df <- z$df

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
  ols <- attr(x, "ols")
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
                 if (ols) signif.legend = FALSE else signif.legend = TRUE,
                 has.Pvalue = TRUE, P.values = TRUE, na.print = "NA", ...)


    if (ols) {
      cat("\nCoefficients OLS:\n")

      printCoefmat(coefs_ols, digits = digits, signif.stars = signif.stars,
                   signif.legend = TRUE, has.Pvalue = TRUE, P.values = TRUE,
                   na.print = "NA", ...)
    }

    cat("\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
        "on", x$df[1], "DF, p-value",
        format.pval(x$pv_waldjoint, digits = digits))
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
  delta <- object[["delta"]]
  if (mat) {
    # upperTriangle(delta, diag = diag) <- NA
  }else{
    delta <- drop(delta)
  }
  delta
}

varcov <- function(object) {
  if (!inherits(object, c("ivx", "summary.ivx"))) stop("Wrong object", call. = F)
  varcov <- drop(object[["varcov"]])
  varcov
}

coefficients.ivx <- function(x, ols = FALSE) {
  if (ols) {
    x$coefficients
  }else{
    x$coefficients_ols
  }
}

#' Fitting IVX-AR Models
#'
#' ivx_ar implements the Yang et al (2020) new instrumental variable based Wald statistic
#' (IVX-AR) which accounts for serial correlation and heteroscedasticity in the error
#' terms of the linear predictive regression model.
#'
#' @inheritParams ivx
#' @param ar Method to include the autoregressive terms. "auto" find the optimal
#' ar order by using the information criteria. \code{ar = 0} reduces to simple \code{\link{ivx}}.
#' \code{ar > 1} uses a fixed order to estimate the model.
#' @param ar_ic Information criterion to be used in model selection.
#' @param ar_max	 Maximum ar order of model to fit.
#' @param ar_grid The ar grid sequence of which to iterate.
#'
#' @references Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the
#' Predictability of US Housing Price Index Returns Based on an IVX-AR Model.
#' Journal of the American Statistical Association, 1-22. DOI:
#' \doi{10.1080/01621459.2019.1686392}
#'
#' @export
#' @examples
#'
#' ivx_ar(hpi ~ log(res) + cpi, ylpc)
#'
#' ivx_ar(hpi ~ log(res) + cpi, ylpc, ar_ic = "aic")
#'
#' ivx_ar(hpi ~ log(res) + cpi, ylpc, ar = 1)
#'
ivx_ar <- function(formula, data, horizon, ar = "auto", ar_ic = c("bic", "aic", "aicc"),
                   ar_max = 5, ar_grid = function(x) seq(x - 0.3, x + 0.3, by = 0.02),
                   na.action, contrasts = NULL, offset, model = TRUE, x = FALSE, y = FALSE,
                   ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  if (missing(horizon)) {
    horizon <- cl$horizon <- 1
  }
  if (!ar %in% c("auto", "forecast", c(0:24))) {
    stop("`ar` should be either 'auto' or a non-negative integer.",
      call. = FALSE
    )
  }
  ar_ic <- match.arg(ar_ic)
  if (ar_max != trunc(ar_max) || ar_max <= 0) {
    stop("`ar_max` should be a positive integer.", call. = FALSE)
  }
  if (!is.function(ar_grid)) {
    stop("`ar_grid` should be function with sequence generation see `?seq`",
      call. = FALSE
    )
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "horizon", "na.action", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf["horizon"] <- NULL
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0) {
    warning("ivx estimation does not include an intercept by construction",
      call. = FALSE
    )
  }
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  if (is.matrix(y)) {
    stop("multivariate model are not available", call. = FALSE)
  }
  ny <- length(y)
  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)
  z <- ivx_ar_fit(y, x,
    horizon = horizon, ar = ar, ar_max = ar_max, ar_ic = ar_ic,
    ar_grid = ar_grid, offset = offset, ...
  )
  class(z) <- if (ar == 0) "ivx" else c("ivx_ar", "ivx")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$assign <- attr(x, "assign")
  if (model) {
    z$model <- mf
  }
  if (ret.x) {
    z$x <- x
  }
  if (ret.y) {
    z$y <- y
  }
  z
}



#' Fitter Functions for IVX-AR Models
#'
#' Basic function called by `ivx_ar` to fit predictive models.
#' These should only be used directly by experienced users.
#'
#' @inheritParams stats::lm.fit
#' @inheritParams ivx_ar
#' @param ... Further arguments passed to the function which is fitting the best AR model.
#' If `ar = "auto"` then the internal function `auto_ar` is used, if `ar = "forecast"` then
#' the the function `forecast::auto.arima` is used. If ar is of fixed length then `arima` is used.
#'
#' @export
#' @importFrom stats arima var
#' @examples
#'
#' ivx_ar_fit(monthly$Ret, as.matrix(monthly$LTY))
#'
#' ivx_ar_fit(monthly$Ret, as.matrix(monthly$LTY), ar = 1)
#'
ivx_ar_fit <- function(y, x, horizon = 1, offset = NULL, ar = "auto", ar_max = 5, ar_ic = "bic",
                       ar_grid = function(x) seq(x - 0.3, x + 0.3, by = 0.02), ...) {
  mdl_ivx <- ivx_fit(y, x, horizon = horizon)
  if (ar == "auto") {
    mdl_ar <- auto_ar(mdl_ivx$ols$residuals, d = 0, max.p = ar_max, ar_ic = ar_ic, ...)
  } else if (ar == "forecast") {
    requireNamespace("forecast", quietly = TRUE)
    mdl_ar <- forecast::auto.arima(mdl_ivx$ols$residuals, d = 0, max.p = ar_max, max.q = 0, ic = ar_ic, ..)
  } else if (ar == 0) {
    message("Using `ivx` instead.")
    return(mdl_ivx)
  } else {
    mdl_ar <- arima2(mdl_ivx$ols$residuals, order = c(ar, 0, 0), include.mean = FALSE, ...)
  }
  ar_coefs <- coefficients(mdl_ar)
  # in case the arima does not converge
  if (is_numeric0(ar_coefs)) {
    return(
      list(
        coefficients = numeric(), residuals = y,
        fitted = 0 * y, df.residuals = length(y),
        ar_method = ar, ar_aic = ar_ic
      )
    )
  }

  res_ar <- residuals(mdl_ar)
  q <- length(ar_coefs)
  grid_seq <- sapply(ar_coefs, ar_grid)
  ngrid <- nrow(grid_seq)

  res_ivx <- vector("list", length = ngrid)
  rse <- vector("numeric", length = ngrid)
  for (i in 1:ngrid) {
    y_adj <- tilt(y, grid_seq[i, ], q)
    x_adj <- tilt(x, grid_seq[i, ], q)
    res_ivx[[i]] <- ivx::ivx_fit(y_adj, x_adj, horizon = horizon)
    eps <- y_adj - sum(x_adj * res_ivx[[i]]$coefficients)
    rse[i] <- var(eps[!is.infinite(eps)])
  }

  rse_min <- which.min(rse)
  z <- res_ivx[[rse_min]]
  z$rse <- rse[rse_min]
  z$coefficients_ar <- grid_seq[rse_min, ]
  z$ar_method <- ar
  z$ar_ic <- ar_ic

  z$Wald_AR <- ac_test_wald(mdl_ivx$ols$residuals, q)
  z$q <- q
  z
}


#' @rdname ivx_ar
#' @param x an object of class "ivx_ar", usually, a result of a call to ivx_ar.
#' @inheritParams stats::summary.lm
#' @export
print.ivx_ar <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    sep = ""
  )
  cat("\n\nLag Selection:\n",
    if (x$ar_method == "auto") paste0("Auto (", x$ar_ic, ")") else "Fixed",
    " with AR terms q = ", x$q, "\n\n",
    sep = ""
  )
  res <- x$coefficients
  if (length(res)) {
    cat("Coefficients:\n")
    print.default(format(res, digits = digits), print.gap = 2L, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}


#' Summarizing IVX-AR Model Fits
#'
#' summary method for class "ivx".
#'
#' @param object  object of class "ivx_ar", usually, a result of a call to ivx_ar.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @importFrom stats printCoefmat
#' @importFrom stats pt qchisq
#' @examples
#' mod <- ivx_ar(Ret ~ LTY, data = kms)
#'
#' summary(mod)
summary.ivx_ar <- function(object, ...) {
  z <- object

  if (is.null(z$terms)) {
    stop("invalid 'ivx' object: no 'terms' components")
  }
  if (!inherits(object, "ivx_ar")) {
    stop("calling summary.ivx(<fake-ivx-object>) ...")
  }

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
  ans$pv_waldjoint <- 1 - pchisq(z$Wald_Joint, z$df)
  ans$df <- z$df

  ans$df.residuals <- z$df.residuals
  rss <- sum(ans$residuals^2)
  mss <- sum(ans$fitted^2)
  n <- NROW(ans$residuals)
  ans$r.squared <- mss / (rss + mss)
  ans$adj.r.squared <- 1 - (1 - ans$r.squared) * n / ans$df.residuals

  ans$q <- z$q
  ans$chi_crit <- qchisq(0.95, ans$q)
  ans$Wald_AR <- z$Wald_AR
  ans$pv_waldar <- 1 - pchisq(ans$Wald_AR, ans$q)

  ans$rse <- z$rse
  ans$ar_coefs <- z$ar_coefs
  ans$ar_method <- z$ar_method
  ans$ar_aic <- z$ar_aic

  if (is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.ivx_ar"

  ans
}

#' @inheritParams stats::summary.lm
#' @rdname summary.ivx_ar
#' @export
print.summary.ivx_ar <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  cat("\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    sep = ""
  )
  cat("\n\n",
    if (x$ar_method == "auto") paste0("Auto (", x$ar_aic, ")") else "Fixed",
    " with AR terms q = ", x$q, "\n\n",
    sep = ""
  )

  if (length(x$aliased) == 0L) {
    cat("No coefficients\n")
  } else {
    coefs_ols <- x$coefficients_ols
    coefs <- x$coefficients
    aliased <- x$aliased

    if (!is.null(aliased) && any(aliased)) {
      cn <- names(aliased)
      civx <- x$coefficients_ols
      coefs <- matrix(NA, NROW(civx), 5, dimnames = list(cn, colnames(civx)))
      coefs[!aliased, ] <- civx
    }

    cat("Coefficients:\n")

    printCoefmat(coefs,
      digits = digits, signif.stars = signif.stars,
      signif.legend = TRUE, has.Pvalue = TRUE, P.values = TRUE,
      na.print = "NA", ...
    )

    cat(
      "\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
      "on", x$df, "DF, p-value",
      format.pval(x$pv_waldjoint, digits = digits)
    )

    cat("\nMultiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits))

    cat(
      "\nWald AR statistic:", formatC(x$Wald_AR, digits = digits),
      "on", x$q, "DF, p-value",
      format.pval(x$pv_waldar, digits = digits)
    )

    cat("\n")
  }
  cat("\n")
  invisible(x)
}

#' @importFrom stats embed
tilt <- function(x, grid_vec, ar_length) {
  x <- as.matrix(x)
  out <- matrix(NA, nrow = NROW(x) - ar_length, ncol = NCOL(x))
  for (i in 1:NCOL(x)) {
    ds <- embed(x[, i], ar_length + 1)
    out[, i] <- ds[, 1] - ds[, -1, drop = FALSE] %*% grid_vec
  }
  colnames(out) <- colnames(x)
  out
}

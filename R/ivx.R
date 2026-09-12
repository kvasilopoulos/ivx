#' Fitting IVX Models
#'
#' ivx fits predictive regression models. The method allows standard
#' chi-square testing for regressors with different degrees of persistence,
#' from stationary to mildly explosive, and can be used for both short-
#' and long-horizon predictive regressions.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data n optional data frame, list or environment (or object coercible by
#' \code{\link[base:as.data.frame]{as.data.frame}} to a data frame) containing
#' the variables in the model. If not found in data, the variables are taken
#' from environment(formula), typically the environment from which lm is called.
#' @param horizon is the horizon (default horizon = 1 corresponds to a
#' short-horizon regression). For `horizon > 1` the estimator and modified
#' Wald statistic are those of Kostakis et al. (2023), eqs (15) and (23):
#' K-period sums of the response and of the lagged predictors, a single-lag
#' instrument in the signal matrix and the K-period sum of the instrument in
#' the covariance.
#' @param weights an optional vector of weights to be used in the fitting process.
#' Should be `NULL` or a numeric vector. If non-NULL, weighted least squares is used
#' with weights `weights` (that is, minimizing `sum(w*e^2)`); otherwise ordinary
#' least squares is used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is set by the na.action setting of \code{\link[base:options]{options}},
#' and is \code{\link[stats:na.fail]{na.fail}} if that is unset. The ‘factory-fresh’
#' default is \code{\link[stats:na.fail]{na.omit}}. Another possible value is \code{NULL},
#' no action. Value \code{\link[stats:na.fail]{na.exclude}} can be useful.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{\link[stats:model.matrix]{model.matrix.default}}.
#' @param offset 	this can be used to specify an a priori known component to be
#' included in the linear predictor during fitting. This should be NULL or a
#' numeric vector or matrix of extents matching those of the response. One or
#' more offset terms can be included in the formula instead or as well, and if more
#' than one are specified their sum is used. See \link[stats:model.extract]{model.offset}
#' @param ... additional arguments to be passed to the low level regression fitting
#' functions (see \link[stats:lm]{lm}).
#' @param model logical. If `TRUE` the model.frame of the fit is returned.
#' @param x logical. If `TRUE` the model.matrix of the fit is returned.
#' @param y logical. If `TRUE` the response of the fit is returned.
#' @param beta,cz tuning parameters of the IVX instrument
#' \eqn{z_t = \sum_{j=0}^{t-1} (1 - c_z/n^eta)^j \Delta x_{t-j}}. Defaults
#' (`beta = 0.95`, `cz = 1`) follow Kostakis et al. (2015).
#' @param bandwidth Newey-West bandwidth for the long-run covariance estimate.
#' The default `NULL` uses \eqn{n^{1/3}} as in Kostakis et al. (2015).
#' @param robust logical. If `TRUE` the Eicker-White (heteroskedasticity-robust)
#' form of the IVX covariance matrix is used (Demetrescu et al., 2023). Only
#' available for `horizon = 1`.
#' @param lag_y logical. If `TRUE` the regression is augmented with the lagged
#' dependent variable (column `y_lag`), instrumented by itself, as in
#' Demetrescu (2014): this can raise the local power of the IVX test when the
#' predictors are highly persistent and endogenous, at no cost otherwise. The
#' joint Wald statistic still tests only the predictors. Only for `horizon = 1`.
#'
#' @return an object of class "ivx".
#'
#' @references Magdalinos, T., & Phillips, P. (2009). Limit Theory for Cointegrated
#' Systems with Moderately Integrated and Moderately Explosive Regressors.
#' Econometric Theory, 25(2), 482-526.
#' @references Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2014).
#' Robust econometric inference for stock return predictability. The Review of
#' Financial Studies, 28(5), 1506-1553.
#' @references Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
#' (2023). Extensions to IVX methods of inference for return predictability.
#' Journal of Econometrics, 237(2), 105271.
#' @references Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023).
#' Taking stock of long-horizon predictability tests: Are factor returns
#' predictable? Journal of Econometrics, 237(2), 105380.
#' @references Demetrescu, M. (2014). Enhancing the local power of IVX-based
#' tests in predictive regressions. Economics Letters, 124(2), 269-273.
#'
#' @aliases ivx
#'
#' @importFrom stats .getXlevels coef coefficients cor lm model.matrix pf
#' model.offset model.response pchisq qnorm residuals symnum is.empty.model
#'
#' @export
#' @examples
#'
#' # Univariate
#' ivx(Ret ~ LTY, data = kms)
#'
#' # Multivariate
#' ivx(Ret ~ LTY + TBL, data = kms)
#'
#' # Longer horizon
#' ivx(Ret ~ LTY + TBL, data = kms, horizon = 4)
#'
#' wt <- runif(nrow(kms))
#' ivx(Ret ~ LTY, data = kms, weights = wt)
#'
#' # lag-augmented IVX (Demetrescu, 2014)
#' ivx(Ret ~ DP, data = kms, lag_y = TRUE)
#'
ivx <- function(formula, data, horizon, na.action, weights,
                contrasts = NULL, offset, model = TRUE, x = FALSE, y = FALSE,
                beta = 0.95, cz = 1, bandwidth = NULL, robust = FALSE,
                lag_y = FALSE, ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  if (missing(horizon))
    horizon <- cl$horizon <- 1

  ## keep only the arguments which should go into the model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "horizon", "na.action", "weights", "offset"),
             names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame) # was as.name("model.frame"), but
  ##    need "stats:: ..." for non-standard evaluation
  mf["horizon"] <- NULL

  mf <- eval.parent(mf)
  # mf <- eval(mf, parent.frame())

  # if (method == "model.frame") return(mf)

  ## 1) allow model.frame to update the terms object before saving it.
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0) {
    warning("ivx estimation does not include an intercept by construction",
            call. = FALSE
    )
  }
  attr(mt, "intercept") <- 0

  y <- model.response(mf, "numeric")
  if(is.matrix(y)) {  # Disable nultivariate model
    stop("multivariate model is not available",call. = FALSE)
  }
  ny <- length(y)
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) {
    stop("'weights' must be a numeric vector")
  }
  offset <- model.offset(mf)
  if (!is.null(offset)) {
    offset <- as.vector(offset)
    if (NROW(offset) != ny)
      stop(
        gettextf("number of offsets is %d, should equal %d (number of observations)", NROW(offset), ny),
        domain = NA)
  }

  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = numeric(), residuals = y, fitted.values = 0 *y,
              weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=0) else ny)
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w)) {
      ivx_fit(y, x, horizon = horizon, offset = offset, beta = beta, cz = cz,
              bandwidth = bandwidth, robust = robust, lag_y = lag_y, ...)
    }else  {
      ivx_wfit(y, x, w, horizon = horizon, offset = offset, beta = beta, cz = cz,
               bandwidth = bandwidth, robust = robust, lag_y = lag_y, ...)
    }
  }
  class(z) <- "ivx"
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
ivx_fit <- function(y, x, horizon = 1, offset = NULL, beta = 0.95, cz = 1,
                    bandwidth = NULL, robust = FALSE, lag_y = FALSE, ...) {
  n <- NROW(x)
  p <- NCOL(x)

  if (is.null(n)) stop("'x' must be a matrix")
  if (n == 0L) stop("0 (non-NA) cases")
  if (p == 0L) {
    return(list(
      coefficients = numeric(), residuals = y,
      fitted = 0 * y, df.residuals = length(y)
    ))
  }

  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1) {
    y <- drop(y)
  }
  if (!is.null(offset)) {
    y <- y - offset
  }
  if (NROW(y) != n) stop("incompatible dimensions")

  chkDots(...)

  x <- lag_y_augment(x, y, lag_y)
  z <- ivx_fit_cpp(y, x, horizon, beta, cz, bandwidth %||% -1L, robust, as.integer(lag_y))
  ivx_out(z, x, horizon, beta, cz, robust, lag_y)
}

# append y as the last regressor: ivx_fit_cpp lags the regressors, so the
# column enters the regression as y_{t-1} (Demetrescu, 2014)
lag_y_augment <- function(x, y, lag_y) {
  if (!lag_y) return(x)
  asgn <- attr(x, "assign")
  x <- cbind(x, y_lag = y)
  attr(x, "assign") <- c(asgn, max(asgn, 0L) + 1L)
  x
}

# common post-processing of the C++ output for ivx_fit and ivx_wfit
ivx_out <- function(z, x, horizon, beta, cz, robust, lag_y = FALSE) {
  p <- NCOL(x)

  cnames <- colnames(x)
  if (is.null(cnames))
    cnames <- paste0("x", 1L:p)

  # stats IVX
  coef <- drop(z$Aivx)
  se <- sqrt(diag(z$varcov))
  wald_ind <- drop(z$wivxind)
  names(coef) <- names(se) <- names(wald_ind) <- cnames

  # stats OlS
  lols <- list(
    coefficients = drop(z$ols$Aols),
    se = drop(z$ols$se),
    tstat = drop(z$ols$tstat),
    residuals = drop(z$ols$residuals)
  )
  names(lols$coefficients) <- names(lols$se) <- names(lols$tstat) <- c("Intercept", cnames)

  # Iniitial
  linitial <- list(
    fitted = drop(z$initial$fitted),
    intercept = drop(z$initial$intercept)
  )

  # data and datam
  colnames(z$data$X) <- cnames
  colnames(z$datam$X) <- cnames

  list(
    coefficients = coef,
    se = se,
    tstat = coef / se,
    intercept = drop(z$intercept),
    fitted = drop(z$fitted),
    residuals = drop(z$residuals),
    Wald_Joint = drop(z$wivx),
    Wald_Ind = wald_ind,
    rank = z$rank,
    horizon = horizon,
    df.residuals = z$df.residuals,
    df = z$df,
    assign = attr(x, "assign"),
    cnames = cnames,
    AR = data.frame(
      Rn = z$Rn,
      Rz = z$Rz,
      row.names = cnames
    ),
    data = z$data,
    datam = z$datam,
    delta = z$delta,
    vcov = z$varcov,
    robust = robust,
    lag_y = lag_y,
    tuning = list(beta = beta, cz = cz, bandwidth = z$bandwidth),
    ols = lols,
    initial = linitial
  )
}


#' @rdname ivx_fit
ivx_wfit <- function(y, x, w, horizon = 1, offset = NULL, beta = 0.95, cz = 1,
                     bandwidth = NULL, robust = FALSE, lag_y = FALSE, ...) {
  n <- nrow(x)
  ny <- NCOL(x)
  if (is.null(n)) {
    stop("'x' must be a matrix")
  }
  if (n == 0) {
    stop("0 (non-NA) cases")
  }
  if (!is.null(offset)) {
    y <- y - offset
  }
  if (NROW(y) != n | length(w) != n) {
    stop("incompatible dimensions")
  }
  if (any(w < 0 | is.na(w))) {
    stop("missing or negative weights not allowed")
  }

  chkDots(...)
  x.asgn <- attr(x, "assign")
  zero.weights <- any(w == 0)
  if (zero.weights) {
    save.r <- y
    save.f <- y
    save.w <- w
    ok <- w != 0
    nok <- !ok
    w <- w[ok]
    x0 <- x[!ok, , drop = FALSE]
    x <- x[ok, , drop = FALSE]
    n <- nrow(x)
    y0 <- if (ny > 1L) {
      y[!ok, , drop = FALSE]
    } else {
      y[!ok]
    }
    y <- if (ny > 1L) {
      y[ok, , drop = FALSE]
    } else {
      y[ok]
    }
  }
  p <- ncol(x)
  if (p == 0) {
    return(list(
      coefficients = numeric(), residuals = y,
      fitted = 0 * y, weights = w, rank = 0L, df.residuals = length(y)
    ))
  }
  if (n == 0) {
    return(list(
      coefficients = rep(NA_real_, p), residuals = y,
      fitted = 0 * y, weights = w, rank = 0L, df.residuals = 0L
    ))
  }
  wts <- sqrt(w)
  x <- lag_y_augment(x, y, lag_y)
  z <- ivx_fit_cpp(y * wts, x * wts, horizon, beta, cz, bandwidth %||% -1L, robust,
                   as.integer(lag_y))
  out <- ivx_out(z, x, horizon, beta, cz, robust, lag_y)

  out$residuals <- out$residuals / wts[-(1:horizon)]
  out$fitted <- y[-(1:horizon)] - out$residuals

  out$weights <- w
  if (zero.weights) {
    coef <- out$coefficients
    coef[is.na(coef)] <- 0
    f0 <- x0 %*% coef
    if (ny > 1) {
      save.r[ok, ] <- out$residuals
      save.r[nok, ] <- y0 - f0
      save.f[ok, ] <- out$fitted
      save.f[nok, ] <- f0
    }
    else {
      save.r[ok] <- out$residuals
      save.r[nok] <- y0 - f0
      save.f[ok] <- out$fitted
      save.f[nok] <- f0
    }
    out$residuals <- save.r
    out$fitted.values <- save.f
    out$weights <- save.w
  }
  if (!is.null(offset)) {
    out$fitted.values <- out$fitted + offset
  }
  out
}


#' @rdname ivx
#' @param x an object of class "ivx", usually, a result of a call to ivx.
#' @inheritParams stats::summary.lm
#' @export
print.ivx <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
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
summary.ivx <- function(object, ...) {
  z <- object

  if (is.null(z$terms)) {
    stop("invalid 'ivx' object: no 'terms' components")
  }
  if (!inherits(object, "ivx")) {
    stop("calling summary.ivx(<fake-ivx-object>) ...")
  }

  ans <- z[c("call", "terms")]

  ans$aliased <- is.na(z$coefficients)

  p_value_ivx <- 1 - pchisq(z$Wald_Ind, 1)
  ans$coefficients <- cbind(z$coefficients, z$se, z$tstat, z$Wald_Ind, p_value_ivx)
  dimnames(ans$coefficients) <- list(
    z$cnames, c("Estimate", "Std. Error", "t value", "Wald Ind", "Pr(> chi)")
  )
  ans$robust <- z$robust

  ans$vcov <- z$vcov
  dimnames(ans$vcov) <- dimnames(ans$coefficients)[c(1, 1)]

  ans$delta <- z$delta
  colnames(ans$delta) <- z$cnames

  ans$residuals <- z$residuals
  ans$fitted <- z$fitted

  ans$horizon <- z$horizon
  ans$Wald_Joint <- z$Wald_Joint
  ans$df <- z$df
  ans$pv_waldjoint <- 1 - pchisq(z$Wald_Joint, z$df)

  ans$df.residuals <- z$df.residuals
  rss <- sum(ans$residuals^2)
  mss <- sum(ans$fitted^2)
  n <- NROW(ans$residuals)
  ans$r.squared <- mss / (rss + mss)
  ans$adj.r.squared <- 1 - (1 - ans$r.squared) * n / ans$df.residuals


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
                              ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
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
      digits = digits, signif.stars = signif.stars, cs.ind = 1:2, tst.ind = 3:4,
      signif.legend = TRUE, has.Pvalue = TRUE, P.values = TRUE,
      na.print = "NA", ...
    )
    if (isTRUE(x$robust)) cat("(Eicker-White standard errors)
")

    cat(
      "\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
      "on", x$df, "DF, p-value",
      format.pval(x$pv_waldjoint, digits = digits)
    )
    cat("\nMultiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits))
    cat("\n")
  }
  cat("\n")
  invisible(x)
}


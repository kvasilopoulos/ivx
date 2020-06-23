#' Fitting IVX-AR Models
#'
#' ivx_ar implements the Yang et al (2020) new instrumental variable based Wald statistic
#' (IVX-AR) which accounts for serial correlation and heteroscedasticity in the error
#' terms of the linear predictive regression model.
#'
#' @inheritParams ivx
#' @param ar Method to include the autoregressive terms. "auto" find the optimal
#' ar order by using the infomration criteria. \code{ar = 0} reduces to simple \code{\link{ivx}}.
#' \code{ar > 1} uses a fixed order to estimate the model.
#' @param ar_ic Information criterion to be used in model selection.
#' @param ar_max	 Maximum ar order of model to fit.
#' @param ar_grid The ar grid sequence of which to iterate
#'
#' @references Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the
#' Predictability of US Housing Price Index Returns Based on an IVX-AR Model.
#' Journal of the American Statistical Association, 1-22. DOI:
#' \href{https://doi.org/10.1080/01621459.2019.1686392}{10.1080/01621459.2019.1686392}
#'
#'
#' @references Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the predictability
#' of US housing price index returns based on an IVX-AR model. Journal of the American
#' Statistical Association, 1-22.
#'
#' @export
ivx_ar <- function(formula, data, horizon, ar = "auto", ar_ic = c("bic", "aic", "aicc"),
                   ar_max = 5, ar_grid =  function(x) seq(x - 0.3, x + 0.3, by = 0.02),
                   na.action, contrasts = NULL, offset, ...) {
  cl <- match.call()
  if (missing(horizon))
    horizon <- cl$horizon <- 1
  if (! ar %in% c("auto", c(0:50))) {
    stop("`ar` should be either 'auto' or a non-negative integer.",
         call. = FALSE)
  }
  ar_ic <- match.arg(ar_ic)
  if(ar_max != trunc(ar_max) || ar_max <= 0) {
    stop("`ar_max` should be a positive integer.", call. = FALSE)
  }
  # if(!is.function(ar_grid)) {
  #   stop("`ar_grid` should be function with sequence generation see `?seq`",
  #        call. = FALSE)
  # }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "horizon",
               "na.action", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf["horizon"] <- NULL
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0) {
    warning("ivx estimation does not include an intercept by construction",
            call. = FALSE)
  }
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)
  z <- ivx_fit_ar(y, x, horizon = horizon, ar = ar, max = ar_max, ic = ar_ic,
                  grid_seq = ar_grid, offset = offset, ...)
  class(z) <- if (ar==0) "ivx" else "ivx_ar"
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z
}
is_numeric0 <- function(x) {
  is.numeric(x) && length(x) == 0
}


ivx_fit_ar <- function(y, x, horizon = 1, offset = NULL, ar = "auto", max = 5, ic = "bic",
                       grid_seq =  function(x) seq(x - 0.3, x + 0.3, by = 0.02), ...) {

  ivx_mdl <- ivx_fit(y, x, horizon = horizon)
  if (ar == "auto") {
    # ar_mdl <- forecast::auto.arima(ivx_mdl$residuals_ols, d = 0, max.q = 0, ic = ic)
    ar_mdl <- auto_ar(ivx_mdl$residuals_ols, d = 0, max.p = max, ic = ic)
  }else if (ar == 0) {
    message("Using `ivx` instead.")
    return(ivx_mdl)
  }else {
    ar_mdl <- arima(ivx_mdl$residuals_ols, order = c(ar, 0, 0), include.mean = FALSE, method = "CSS")
  }
  ar_coefs <- coefficients(ar_mdl)
  if (is_numeric0(ar_coefs)) {
    return(ar_mdl)
  }

  ar_res <- residuals(ar_mdl)
  q <- length(ar_coefs)
  ar_grid <- sapply(ar_coefs, grid_seq)
  ngrid <- nrow(ar_grid)

  ivx_res <- vector("list", length = ngrid)
  rse <- vector("numeric", length = ngrid)
  for (i in 1:ngrid) {
    y_adj <- tilt(y, ar_grid[i,], q)
    x_adj <- tilt(x, ar_grid[i,], q)
    ivx_res[[i]] <- ivx:::ivx_fit(y_adj, x_adj, horizon = horizon)
    eps <- y_adj - sum(x_adj * ivx_res[[i]]$coefficients)
    rse[i] <- var(eps[!is.infinite(eps)])
  }

  rse_min <- which.min(rse)
  z <- ivx_res[[rse_min]]
  z$rse <- rse[rse_min]
  z$ar_coefs <- ar_grid[rse_min, ]
  z$ar_method <- ar
  z$ar_aic <- ic

  nr <- NROW(ivx_mdl$residuals_ols)
  u <- matrix(ivx_mdl$residuals_ols[(q + 2 - 1):(q + 2 - q)], nrow = q)
  sum_u <- tcrossprod(u)
  sum_v <- tcrossprod(u) * ar_res[q + 2]
  for (t in c(q + 3):nr) {
    uu <- matrix(ivx_mdl$residuals_ols[(t - 1):(t - q)], nrow = q)
    sum_u <- sum_u + tcrossprod(uu)
    sum_v <- sum_v + tcrossprod(uu) *  ar_res[t]^2
  }
  W <- matrix(ar_coefs, nrow = 1) %*% sum_u %*% solve(sum_v) %*% sum_u %*% matrix(ar_coefs, ncol = 1)

  z$Wald_AR <- W
  z$q <- q
  z
}


#' @rdname ivx_ar
#' @param x an object of class "ivx_ar", usually, a result of a call to ivx_ar.
#' @inheritParams stats::summary.lm
#' @export
print.ivx_ar <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\n",
      if(x$ar_method == "auto") paste0("Auto (", x$ar_aic, ")") else "Fixed",
      " with AR terms q = ", x$q, "\n\n", sep ="")
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
#' @importFrom stats pt
#' @examples
#' mod <- ivx_ar(Ret ~ LTY, data = kms)
#'
#' summary(mod)
summary.ivx_ar <- function(object,  ...) {
  z <- object

  if (is.null(z$terms))
    stop("invalid 'ivx' object: no 'terms' components")
  if (!inherits(object, "ivx_ar"))
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
                              ...){

  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\n",
      if(x$ar_method == "auto") paste0("Auto (", x$ar_aic, ")") else "Fixed",
      " with AR terms q = ", x$q, "\n\n", sep ="")

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
    cat("\nWald AR statistic:", formatC(x$Wald_AR, digits = digits),
        "on", x$q, "DF, p-value",
        format.pval(x$pv_waldar, digits = digits))
  }
  cat("\n")
  invisible(x)
}




auto_ar <- function(x, max.p = 5, ic = "aic", ...) {
  p <- max.p
  bestfit <- arima2(x, order = c(p, 0, 0), include.mean = FALSE)
  while (p >= 0) {
    fit <- arima2(x, order = c(p, 0, 0), include.mean = FALSE)
    # print(fit[[ic]])
    if (fit[[ic]] < bestfit[[ic]]) {
      bestfit <- fit
    }
    p <- p - 1
  }
  bestfit
}

arima2 <- function(x, order = c(0, 0, 0), ic = "aic", method = "CSS", ...) {
  fit <- arima(x = x, order = order, ...)
  npar <- length(fit$coef[fit$mask]) + 1
  n <- NROW(x)
  m <- frequency(x)
  nstar <- n - m
  dots <- list(...)
  if (is.null(dots$xreg)) {
    nxreg <- 0
  } else {
    nxreg <- ncol(as.matrix(xreg))
  }
  # if (method == "CSS") {
  #   fit$aic <- nstar * log(fit$sigma2) + 2 * npar
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



tilt <- function(x, grid_vec, ar_length) {
  x <- as.matrix(x)
  out <- matrix(NA, nrow = NROW(x) - ar_length, ncol = NCOL(x))
  for (i in 1:NCOL(x)) {
    ds <- embed(x[,i], ar_length + 1)
    out[,i] <- ds[,1] - ds[,-1, drop = FALSE] %*% grid_vec
  }
  colnames(out) <- colnames(x)
  out
}

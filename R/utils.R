is_numeric0 <- function(x) {
  is.numeric(x) && length(x) == 0
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Model frame, terms, response and design matrix shared by the formula
# interfaces. `call` is the caller's match.call(expand.dots = FALSE) and `env`
# the frame to evaluate it in; `extra` names further model-frame arguments
# (weights, offset). The intercept is always removed from the terms: ivx-type
# estimators demean by construction.
ivx_frame <- function(call, env, contrasts = NULL, extra = character(), univariate = TRUE) {
  mf <- call
  m <- match(c("formula", "data", "na.action", extra), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(mf, env)
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0) {
    warning("ivx estimation does not include an intercept by construction", call. = FALSE)
  }
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  if (univariate && is.matrix(y)) stop("multivariate model is not available", call. = FALSE)
  list(mf = mf, terms = mt, y = y, x = model.matrix(mt, mf, contrasts))
}

# attach the usual model components to a fitted object
ivx_finish <- function(z, fr, call, model = TRUE, ret.x = FALSE, ret.y = FALSE) {
  z$na.action <- attr(fr$mf, "na.action")
  z$contrasts <- attr(fr$x, "contrasts")
  z$xlevels <- .getXlevels(fr$terms, fr$mf)
  z$call <- call
  z$terms <- fr$terms
  if (model) z$model <- fr$mf
  if (ret.x) z$x <- fr$x
  if (ret.y) z$y <- fr$y
  z
}

# ADF regression dx_t = mu + theta x_{t-1} + sum_{i=1}^k psi_i dx_{t-i} + e_t with the
# number of lagged differences k in 0..kmax chosen by BIC or by the MBIC of Ng & Perron
# (2001) on the common sample t = kmax + 2..n, then refitted on t = k + 2..n. Returns
# k, theta, psi, the residuals and the normalised-bias statistic n theta / (1 - sum psi).
adf_lag <- function(x, kmax, ic = c("bic", "mbic")) {
  ic <- match.arg(ic)
  n <- length(x)
  dx <- diff(x)
  fit <- function(k, start) {
    idx <- start:(n - 1)
    X <- cbind(1, x[idx])
    if (k > 0) X <- cbind(X, sapply(seq_len(k), function(i) dx[idx - i]))
    f <- lm.fit(X, dx[idx])
    list(coef = f$coefficients, resid = f$residuals, xl = x[idx])
  }
  crit <- sapply(0:kmax, function(k) {
    f <- fit(k, kmax + 1)
    Te <- length(f$resid)
    s2 <- sum(f$resid^2) / Te
    if (ic == "bic") Te * log(s2) + (k + 2) * log(Te)
    else log(s2) + log(Te) * (f$coef[2]^2 * sum(f$xl^2) / s2 + k) / Te
  })
  k <- which.min(crit) - 1
  f <- fit(k, k + 1)
  psi <- if (k > 0) unname(f$coef[-(1:2)]) else numeric()
  list(k = k, theta = unname(f$coef[2]), psi = psi, resid = f$resid,
       stat = unname(n * f$coef[2] / (1 - sum(psi))))
}

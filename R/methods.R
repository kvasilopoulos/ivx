# New methods -------------------------------------------------------------

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
  if (!inherits(object, c("ivx", "summary.ivx", "ivx_ar", "summary.ivx_ar"))) {
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

#' @export
model.matrix.ivx <- function(object, ...) {
  if (n_match <- match("x", names(object), 0L)) {
    object[[n_match]]
  } else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    if (exists(".GenericCallEnv", inherits = FALSE)) {
      NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    } else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      do.call(
        "model.matrix.default",
        c(
          list(object = object, data = data, contrasts.arg = object$contrasts),
          dots
        )
      )
    }
  }
}

#' @export
model.frame.ivx <- function(formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(
      c("formula", "data", "subset", "weights", "na.action", "offset"),
      names(fcall),
      0L
    )
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) {
      env <- parent.frame()
    }
    eval(fcall, env)
  } else {
    formula$model
  }
}

#' @export
logLik.ivx <- function(object, ...) {
  res <- object$residuals
  p <- object$rank
  w <- object$weights
  if (is.null(w)) {
    w <- rep.int(1, length(res))
  } else if (length(w) != length(res)) {
    w <- w[-(1:(length(w) - length(res)))] # the fit drops the first `horizon` rows
  }
  # zero-weight rows and the lag-dropped rows (NA residuals) do not enter
  keep <- w != 0 & !is.na(res)
  res <- res[keep]
  w <- w[keep]
  N <- length(res)
  # Gaussian log-likelihood at the IVX residuals (no REML: there is no QR of an IV fit)
  val <- 0.5 *
    (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  attr(val, "nall") <- N
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}

#' @export
nobs.ivx <- function(object, ...) {
  attr(logLik(object), "nobs")
}

#' @export
#' @importFrom stats weighted.residuals
deviance.ivx <- function(object, ...) {
  sum(weighted.residuals(object)^2, na.rm = TRUE)
}

#' @export
extractAIC.ivx <- function(fit, scale = 0, k = 2, ...) {
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  RSS <- deviance.ivx(fit)
  dev <- if (scale > 0) {
    RSS / scale - n
  } else {
    n * log(RSS / n)
  }
  c(edf, dev + k * edf)
}

#' @export
#' @importFrom stats case.names weights
case.names.ivx <- function(object, full = FALSE, ...) {
  w <- weights(object)
  dn <- names(residuals(object)) %||% as.character(seq_along(residuals(object)))
  if (full || is.null(w)) {
    dn
  } else {
    dn[w != 0]
  }
}

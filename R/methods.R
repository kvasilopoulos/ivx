
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
model.matrix.ivx <- function (object, ...) {
  if (n_match <- match("x", names(object), 0L))
    object[[n_match]]
  else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    if (exists(".GenericCallEnv", inherits = FALSE))
      NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      do.call("model.matrix.default", c(list(object = object,
                                             data = data, contrasts.arg = object$contrasts),
                                        dots))
    }
  }
}

#' @export
model.frame.ivx <- function (formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"),
                      names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset",
                 "weights", "na.action", "offset"),
               names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env))
      env <- parent.frame()
    eval(fcall, env)
  }
  else formula$model
}

model_matrix <- function(object, raw = FALSE) {
  data <- object$data
  ret <- if(raw) data$X else data$Xm
  rownames(ret) <- 1:NROW(ret)
  ret
}

#' @importFrom stats model.frame
model_frame <- function(object, raw = FALSE) {
  data <- object$data
  ret <- if(raw) cbind(data$y, data$X) else cbind(data$ym, data$Xm)
  colnames(ret) <- colnames(model.frame(object))
  rownames(ret) <- 1:NROW(ret)
  ret
}

#' @export
logLik.ivx <- function (object, REML = FALSE, ...) {

  res <- object$residuals
  p <- object$rank
  N <- length(res)
  w <- object$weights
  if (is.null(w)) {
    w <- rep.int(1, N)
  }
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  N0 <- N
  if (REML)
    N <- N - p
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  if (REML)
    val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
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
  dev <- if (scale > 0)
    RSS/scale - n
  else n * log(RSS/n)
  c(edf, dev + k * edf)
}

#' @export
#' @importFrom stats case.names weights
case.names.ivx <- function (object, full = FALSE, ...) {
  w <- weights(object)
  dn <- names(residuals(object))
  if (full || is.null(w))
    dn
  else dn[w != 0]
}


# Unfinished --------------------------------------------------------------

# TODO predict.ivx <- function() {}
# TODO simulate.ivx <- function() {}

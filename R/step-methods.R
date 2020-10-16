
#' @export
deviance.ivx <- function(object, ...) {
  sum(residuals(object)^2, na.rm = TRUE)
}

#' @export
extractAIC.ivx <- function(fit, scale = 0, k = 2, ..) {
  n <- length(fit$residuals)
  edf <- fit$df[1]
  RSS <- deviance.ivx(fit)
  dev <- if (scale > 0)
    RSS/scale - n
  else n * log(RSS/n)
  c(edf, dev + k * edf)
}

check_exact <- function (object) {
  mss <- sum(object$fitted^2)
  rss <- sum(object$residuals^2)
  if (rss < 1e-10 * mss)
    warning("attempting model selection on an essentially perfect fit is nonsense",
            call. = FALSE)
}

#' @export
model.matrix.ivx <- function (object, ...) {
  if (n_match <- match("x", names(object), 0L)) {
    object[[n_match]]
  } else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    if (exists(".GenericCallEnv", inherits = FALSE))
      NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      do.call("model.matrix.default",
              c(list(object = object, data = data, contrasts.arg = object$contrasts), dots))
    }
  }
}

#' @export
model.frame.ivx <- function (formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "na.action", "offset"), names(fcall), 0L)
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


#' object <- ivx(hpi ~ inc + cpi, data = ylpc)
#' object2 <- lm(hpi ~ inc + cpi, data = ylpc)
#' drop1(object)
#' @export
drop1.ivx <-  function (object, scope, scale = 0, all.cols = TRUE, test = c("none", "Chisq", "F"), k = 2, ...) {
  check_exact(object)
  x <- model.matrix(object)
  offset <- model.offset(model.frame(object))
  n <- nrow(x)
  asgn <- attr(x, "assign")
  tl <- attr(object$terms, "term.labels")
  if (missing(scope)) {
    scope <- drop.scope(object)
  }else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ndrop <- match(scope, tl)
  ns <- length(scope)
  rdf <- object$df.residual
  chisq <- deviance.ivx(object)
  dfs <- numeric(ns)
  RSS <- numeric(ns)
  y <- object$residuals + object$fitted
  na.coef <- seq_along(object$coefficients)[!is.na(object$coefficients)]
  for (i in seq_len(ns)) {
    ii <- seq_along(asgn)[asgn == ndrop[i]]
    jj <- setdiff(if (all.cols) seq(ncol(x)) else na.coef, ii)
    z <- ivx_fit(y, x[-1, jj, drop = FALSE], offset = offset)
    # dfs[i] <- z$rank
    oldClass(z) <- "ivx"
    RSS[i] <- deviance(z)
  }
  # scope <- c("<none>", scope)
  # dfs <- c(object$rank, dfs)
  # RSS <- c(chisq, RSS)
  if (scale > 0) {
    aic <- RSS/scale - n + k * dfs
  } else {
    aic <- n * log(RSS/n) + k * dfs
  }
  dfs <- dfs[1L] - dfs
  dfs[1L] <- NA
  aod <- data.frame(
    Df = dfs,
    `Sum of Sq` = c(NA, RSS[-1L] - RSS[1L]),
    RSS = RSS, AIC = aic,
    row.names = scope, check.names = FALSE)
  if (scale > 0)
    names(aod) <- c("Df", "Sum of Sq", "RSS","Cp")
  test <- match.arg(test)
  # if (test == "Chisq") {
  #   dev <- aod$"Sum of Sq"
  #   if (scale == 0) {
  #     dev <- n * log(RSS/n)
  #     dev <- dev - dev[1L]
  #     dev[1L] <- NA
  #   }
  #   else dev <- dev/scale
  #   df <- aod$Df
  #   nas <- !is.na(df)
  #   dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail = FALSE)
  #   aod[, "Pr(>Chi)"] <- dev
  # }
  # else if (test == "F") {
  #   dev <- aod$"Sum of Sq"
  #   dfs <- aod$Df
  #   rdf <- object$df.residual
  #   rms <- aod$RSS[1L]/rdf
  #   Fs <- (dev/dfs)/rms
  #   Fs[dfs < 1e-04] <- NA
  #   P <- Fs
  #   nas <- !is.na(Fs)
  #   P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
  #   aod[, c("F value", "Pr(>F)")] <- list(Fs,
  #                                         P)
  # }
  head <- c("Single term deletions", "\nModel:",
            deparse(formula(object)), if (scale > 0) paste("\nscale: ",
                                                           format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod

}

#' @export
add1.ivx <- function (object, scope, scale = 0, test = c("none", "Chisq", "F"), x = NULL, k = 2, ...) {
  Fstat <- function(table, RSS, rdf) {
    dev <- table$"Sum of Sq"
    df <- table$Df
    rms <- (RSS - dev)/(rdf - df)
    Fs <- (dev/df)/rms
    Fs[df < .Machine$double.eps] <- NA
    P <- Fs
    nnas <- !is.na(Fs)
    P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas],
                       lower.tail = FALSE)
    list(Fs = Fs, P = P)
  }
  check_exact(object)
  if (missing(scope) || is.null(scope))
    stop("no terms in scope")
  if (!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope))
    stop("no terms in scope for adding to object")
  oTerms <- attr(object$terms, "term.labels")
  int <- attr(object$terms, "intercept")
  y <- object$residuals + object$fitted
  ns <- length(scope)
  dfs <- numeric(ns + 1L)
  names(dfs) <- c("<none>", scope)
  RSS <- dfs
  add.rhs <- eval(str2lang(paste("~ . +", paste(scope, collapse = "+"))))
  new.form <- update.formula(object, add.rhs)
  Terms <- terms(new.form)
  if (is.null(x)) {
    fc <- object$call
    fc$formula <- Terms
    fob <- list(call = fc, terms = Terms)
    class(fob) <- oldClass(object)
    m <- model.frame(fob, xlev = object$xlevels)
    x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- model.offset(m)
    # wt <- model.weights(m)
    oldn <- length(y)
    y <- model.response(m, "numeric")
    newn <- length(y)
    if (newn < oldn)
      warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit",
                               "using the %d/%d rows from a combined fit"),
                      newn, oldn), domain = NA)
  }
  else {
    # wt <- object$weights
    offset <- object$offset
  }
  n <- nrow(x)
  Terms <- attr(Terms, "term.labels")
  asgn <- attr(x, "assign")
  ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
  if (int)
    # ousex[1L] <- TRUE
  iswt <- !is.null(wt)
  X <- x[, ousex, drop = FALSE]
  z <- ivx_fit(y, X, offset = offset)
  dfs[1L] <- z$rank
  class(z) <- "ivx"
  RSS[1L] <- deviance(z)
  sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE),
                   function(x) paste(sort(x), collapse = ":"))
  for (tt in scope) {
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
    usex <- match(asgn, match(stt, sTerms), 0L) > 0L
    X <- x[, usex | ousex, drop = FALSE]
    z <- if (iswt)
      lm.wfit(X, y, wt, offset = offset)
    else lm.fit(X, y, offset = offset)
    class(z) <- "ivx"
    dfs[tt] <- z$rank
    RSS[tt] <- deviance(z)
  }
  if (scale > 0)
    aic <- RSS/scale - n + k * dfs
  else aic <- n * log(RSS/n) + k * dfs
  dfs <- dfs - dfs[1L]
  dfs[1L] <- NA
  aod <- data.frame(
    Df = dfs,
    `Sum of Sq` = c(NA, RSS[1L] -RSS[-1L]),
    RSS = RSS,
    AIC = aic,
    row.names = names(dfs),
    check.names = FALSE)
  if (scale > 0)
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  # if (test == "Chisq") {
  #   dev <- aod$"Sum of Sq"
  #   if (scale == 0) {
  #     dev <- n * log(RSS/n)
  #     dev <- dev[1L] - dev
  #     dev[1L] <- NA
  #   }
  #   else dev <- dev/scale
  #   df <- aod$Df
  #   nas <- !is.na(df)
  #   dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail = FALSE)
  #   aod[, "Pr(>Chi)"] <- dev
  # }
  # else if (test == "F") {
  #   rdf <- object$df.residual
  #   aod[, c("F value", "Pr(>F)")] <- Fstat(aod,
  #                                          aod$RSS[1L], rdf)
  # }
  head <- c("Single term additions", "\nModel:", deparse(formula(object)),
            if (scale > 0) paste("\nscale: ",format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}


#' @export
df.residual.ivx <- function(object, ...) {
  object$df[2]
}


#' @export
anova.ivx <- function (object, ...) {
  if (length(list(object, ...)) > 1L)
    return(anova.lmlist(object, ...))
  if (!inherits(object, "ivx"))
    warning("calling anova.ivx(<fake-ivx-object>) ...")
  ssr <- sum(object$residuals^2)
  mss <- sum(object$fitted^2)
  if (ssr < 1e-10 * mss)
    warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  dfr <- df.residual(object)
  p <- object$rank
  if (p > 0L) {
    p1 <- 1L:p
    comp <- object$effects[p1]
    asgn <- object$assign[qr.lm(object)$pivot][p1]
    nmeffects <- c("(Intercept)", attr(object$terms,
                                       "term.labels"))
    tlabels <- nmeffects[1 + unique(asgn)]
    ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
    df <- c(lengths(split(asgn, asgn)), dfr)
  }
  else {
    ss <- ssr
    df <- dfr
    tlabels <- character()
  }
  ms <- ss/df
  f <- ms/(ssr/dfr)
  P <- pf(f, df, dfr, lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f, P)
  table[length(P), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Residuals"),
                          c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  if (attr(object$terms, "intercept"))
    table <- table[-1, ]
  structure(table, heading = c("Analysis of Variance Table\n",
                               paste("Response:", deparse(formula(object)[[2L]]))),
            class = c("anova", "data.frame"))
}

#' @export
logLik.ivx <- function (object, REML = FALSE, ...) {

  res <- object$residuals
  p <- object$rank
  N <- length(res)
  w <- rep.int(1, N)
  N0 <- N
  if (REML){
    N <- N - p
  }
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  if (REML){
    val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
  }
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}

# methods(class = "ivx")


# y<-rnorm(100)
# x<-rnorm(100)
# m<-lm(y ~ x)
#
# res <- m$residuals
# n <- nrow(m$model)
# w <- rep(1,n) #not applicable
#
# ll < -0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
# ll - logLik(m)==0 #TRUE
#
#
# k.original < -length(m$coefficients)
# df.ll < -k.original+1
#
# bic <- -2 * ll + log(n) * df.ll
# bic - BIC(m)==0 #TRUE
#
# aic <- -2 * ll + 2 * df.ll
# aic - AIC(m)==0 #TRUE
#
#
# extractAIC(lm(y ~ x)) - extractAIC(lm(y ~ 1))
# AIC(lm(y ~ x)) - AIC(lm(y ~ 1))

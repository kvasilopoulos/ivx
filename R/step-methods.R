# Bayesian Criterion k = log(nrow(model.matrix(object)))

assert_ivx <- function(object) {
  if(length(class(object)) > 1) {
    stop("only supported for ivx models")
  }
}


#' @export
#' @importFrom stats model.frame drop.scope terms update.formula deviance formula
drop1.ivx <- function (object, scope, scale = 0, all.cols = TRUE, test = c("none", "Chisq", "F"), k = 2, ...) {

  assert_ivx(object)
  check_exact(object)
  x <- model.matrix(object)
  offset <- model.offset(model.frame(object))
  wt <- object$weights
  iswt <- !is.null(wt)
  n <- nrow(x)
  asgn <- attr(x, "assign")
  tl <- attr(object$terms, "term.labels")
  if (missing(scope)) {
    scope <- drop.scope(object)
  }
  else {
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
    sub_i <- 1 + if(!is.null(object$q)) object$q else 0
    z <- if (iswt) {
      ivx_wfit(y, x[-(1:sub_i), jj, drop = FALSE], wt, offset = offset)
    } else {
      ivx_fit(y, x[-(1:sub_i), jj, drop = FALSE],  offset = offset)
    }
    dfs[i] <- z$rank #if(is.null(z$rank)) 1 else z
    oldClass(z) <- "ivx"
    RSS[i] <- deviance(z)
  }
  scope <- c("<none>", scope)
  dfs <- c(object$rank, dfs)
  RSS <- c(chisq, RSS)
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
    RSS = RSS,
    AIC = aic,
    row.names = scope,
    check.names = FALSE)
  if (scale > 0)
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- aod$"Sum of Sq"
    if (scale == 0) {
      dev <- n * log(RSS/n)
      dev <- dev - dev[1L]
      dev[1L] <- NA
    } else {
      dev <- dev/scale
    }
    df <- aod$Df
    nas <- !is.na(df)
    dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }else if (test == "F") {
    dev <- aod$"Sum of Sq"
    dfs <- aod$Df
    rdf <- object$df.residual
    rms <- aod$RSS[1L]/rdf
    Fs <- (dev/dfs)/rms
    Fs[dfs < 0.0001] <- NA
    P <- Fs
    nas <- !is.na(Fs)
    P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
    aod[, c("F value", "Pr(>F)")] <- list(Fs, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

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

#' @export
#' @importFrom stats arima add.scope update.formula terms model.frame
#' model.weights deviance formula
add1.ivx <- function (object, scope, scale = 0, test = c("none", "Chisq", "F"), x = NULL, k = 2, ...) {
  assert_ivx(object)
  check_exact(object)
  if (missing(scope) || is.null(scope))
    stop("no terms in scope", call. = FALSE)
  if (!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope))
    stop("no terms in scope for adding to object", call. = FALSE)
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
    wt <- model.weights(m)
    oldn <- length(y)
    y <- model.response(m, "numeric")
    newn <- length(y)
    if (newn < oldn)
      warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit",
                               "using the %d/%d rows from a combined fit"),
                      newn, oldn), domain = NA)
  }
  else {
    wt <- object$weights
    offset <- object$offset
  }
  n <- nrow(x)
  Terms <- attr(Terms, "term.labels")
  asgn <- attr(x, "assign")
  ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
  if (int)
    ousex[1L] <- TRUE
  iswt <- !is.null(wt)
  X <- x[, ousex, drop = FALSE]
  z <- if (iswt)
    ivx_wfit(y, X, wt, offset = offset)
  else ivx_fit(y, X, offset = offset)
  dfs[1L] <- z$rank
  class(z) <- "ivx"
  RSS[1L] <- deviance(z)
  sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x), collapse = ":"))
  for (tt in scope) {
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
    usex <- match(asgn, match(stt, sTerms), 0L) > 0L
    X <- x[, usex | ousex, drop = FALSE]
    z <- if (iswt)
      ivx_wfit(y, X, wt, offset = offset)
    else ivx_fit(y, X, offset = offset)
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
    `Sum of Sq` = c(NA, RSS[1L] - RSS[-1L]),
    RSS = RSS, AIC = aic,
    row.names = names(dfs),
    check.names = FALSE)
  if (scale > 0)
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- aod$"Sum of Sq"
    if (scale == 0) {
      dev <- n * log(RSS/n)
      dev <- dev[1L] - dev
      dev[1L] <- NA
    }
    else dev <- dev/scale
    df <- aod$Df
    nas <- !is.na(df)
    dev[nas] <- safe_pchisq(dev[nas], df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }
  else if (test == "F") {
    rdf <- object$df.residual
    aod[, c("F value", "Pr(>F)")] <- Fstat(aod, aod$RSS[1L], rdf)
  }
  head <- c("Single term additions", "\nModel:", deparse(formula(object)),
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
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


# non-exporeted functions from stats ----------------------------------------


check_exact <- function (object) {
  w <- object$weights
  if (is.null(w)) {
    mss <- sum(object$fitted.values^2)
    rss <- sum(object$residuals^2)
  }
  else {
    mss <- sum(w * object$fitted.values^2)
    rss <- sum(w * object$residuals^2)
  }
  if (rss < 0.0000000001 * mss)
    warning("attempting model selection on an essentially perfect fit is nonsense",
            call. = FALSE)
}

safe_pchisq <- function (q, df, ...) {
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

safe_pf <- function (q, df1, ...) {
  df1[df1 <= 0] <- NA
  pf(q = q, df1 = df1, ...)
}



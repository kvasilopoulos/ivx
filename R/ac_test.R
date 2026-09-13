
# AR Wald Test -------------------------------------------------------------

#' Tests  for autocorrelation
#'
#' @description
#'
#' * `ac_test`: For all tests
#' * `ac_test_wald`: Wald test
#' * `ac_test_lb`: Ljung-Box
#' * `ac_test_bp`:  Box-Pierce
#' * `ac_test_bg`: Breusch-Godfrey
#'
#' @param x an `ivx` model or a `numeric vector`, usually the residuals from an ols regression.
#' @param lag the number of lags.
#'
#'
#' @details If p-value < 0.051: You can reject the null hypothesis assuming a
#' 5% chance of making a mistake. So you can assume that your values are showing
#' dependence on each other.
#'
#' @seealso `Box.test` `lmtest::bgtest`
#'
#' @return a numeric scalar or numeric vector.
#'
#' @importFrom stats arima Box.test model.frame
#' @examples
#'
#' mdl <- ivx(hpi ~ cpi + inv, data = ylpc)
#' ac_test_wald(mdl)
#'
#' ac_test(mdl)
#'
#' @name ac_test_
#' @export
ac_test_wald <- function(x, lag) {
  UseMethod("ac_test_wald")
}

#' @export
ac_test_wald.default <- function(x, lag = 1) {
  res_ols <- x
  nr <- NROW(res_ols)
  max_lag <- max(lag)
  war <- vector("numeric", max_lag)
  for(q in 1:max_lag) {
    fit_ar <- arima(res_ols, order = c(q, 0, 0), method = "ML", include.mean = F)
    res_ar <- as.matrix(residuals(fit_ar))
    coef_ar <- coefficients(fit_ar)

    if(q == 1) {
      xx <- coef_ar * sum(res_ols^2) / (nr - 2)
      ss <- sum(res_ols^2 * res_ar^2) / (nr - 2)
      war[1] <- (nr - 2) * xx / ss * xx
    }else{
      u1 <- matrix(res_ols[(q + 2 - 1):(q + 2 - q)], nrow = q)
      sum_u <- tcrossprod(u1)
      sum_v <- tcrossprod(u1) * res_ar[q + 2]

      for (t in c(q + 3):nr) {
        u <- matrix(res_ols[(t - 1):(t - q)], nrow = q)
        sum_u <- sum_u + tcrossprod(u)
        sum_v <- sum_v + tcrossprod(u) * res_ar[t]^2
      }
      war[q] <- matrix(coef_ar, nrow = 1) %*% sum_u %*% solve(sum_v) %*% sum_u %*% matrix(coef_ar, ncol = 1)
    }
  }
  out <- war[lag]
  names(out) <- lag
  class(out) <- "ac_test_"
  attr(out, "pval") <- 1 - pchisq(out, lag)
  attr(out, "tname") <- "Wald"
  out
}

# the tests read the OLS residuals that only ivx() stores; the extensions
# (ivx_ar, ivx_ra, ivx_qr, ivx_iv, arm) inherit the class but not the field
ols_residuals <- function(x) {
  if (class(x)[1L] != "ivx") {
    stop("ac_test is defined for ivx() fits only, not ", class(x)[1L], call. = FALSE)
  }
  x$ols$residuals
}

#' @export
ac_test_wald.ivx <- function(x, lag = 1) {
  res <- ols_residuals(x)
  ac_test_wald.default(res, lag = lag)
}

#' @export
print.ac_test_ <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  tname <- attr(x, "tname")
  pval <- attr(x, "pval")

  stars <- stars_pval(pval)
  stats <- paste0(formatC(x, digits = digits), stars)

  out <- data.frame(names(x),stats)
  colnames(out) <- c("Lag", tname)
  print(out, row.names = FALSE)
}


# Ljung-Box ---------------------------------------------------------------

#' @rdname ac_test_
#' @param lag the number of lags.
#' @export
ac_test_lb <- function(x, lag) {
  UseMethod("ac_test_lb")
}

#' @export
ac_test_lb.default <- function(x, lag = 1) {
  max_lag <- max(lag)
  lb <- vector("numeric", max_lag)
  for(i in 1:max_lag) {
    lb[i] <- Box.test(x, lag = i, type = "Ljung-Box")$statistic
  }
  out <- lb[lag]
  names(out) <- lag
  attr(out, "pval") <- 1 - pchisq(out, lag)
  class(out) <- "ac_test_"
  attr(out, "tname") <- "Ljung-Box"
  out
}

#' @export
ac_test_lb.ivx <- function(x, lag = 1) {
  res <- ols_residuals(x)
  ac_test_lb.default(res, lag)
}

# Box-Pierce --------------------------------------------------------------

#' @rdname ac_test_
#' @export
ac_test_bp <- function(x, lag) {
  UseMethod("ac_test_bp")
}

#' @export
ac_test_bp.default <- function(x, lag = 1) {
  max_lag <- max(lag)
  bp <- vector("numeric", max_lag)
  for(i in 1:max_lag) {
    bp[i] <- Box.test(x, lag = i, type = "Box-Pierce")$statistic
  }
  out <- bp[lag]
  names(out) <- lag
  attr(out, "pval") <- 1 - pchisq(out, lag)
  class(out) <- "ac_test_"
  attr(out, "tname") <- "Box-Pierce"
  out
}

#' @export
ac_test_bp.ivx <- function(x, lag = 1) {
  res <- ols_residuals(x)
  ac_test_bp.default(res, lag)
}

# Breusch–Godfrey ----------------------------------------------------------

#' @param order lag TODO
#' @param type the type of test statistic to be returned. Either "Chisq" for
#' the Chi-squared test statistic or "F" for the F test statistic.
#' @param fill starting values for the lagged residuals in the auxiliary regression.
#' By default 0 but can also be set to NA.
#'
#' @rdname ac_test_
#' @export
ac_test_bg <- function(x, order, type, fill) {
  UseMethod("ac_test_bg")
}


#' @export
ac_test_bg.ivx <- function(x, order = 1, type = c("Chisq", "F"), fill = 0) {

  X <- model.frame(x)
  Y <- model.matrix(x)
  order <- 1:order
  n <- nrow(X)
  k <- ncol(X)
  m <- length(order)

  res_ <- c(ols_residuals(x), rep(0, x$horizon))

  Z <- sapply(order, function(x) c(rep(fill, length.out = x), res_[1:(n - x)]))
  auxfit <- lm(res_ ~.,  cbind(res_, X, Z))

  type <- match.arg(type)
  switch(
    type,
    Chisq = {
      bg <- n * sum(auxfit$fitted.values^2)/sum(res_^2)
      p.val <- 1 - pchisq(bg, m)
      df <- m
      names(df) <- "df"
    }, F = {
      uresi  <- auxfit$residuals
      bg <- ((sum(res_^2) - sum(uresi^2))/m)/(sum(uresi^2)/(n - k - m))
      df <- c(m, n - k - m)
      names(df) <- c("df1", "df2")
      p.val <- 1 - pf(bg, df1 = df[1], df2 = df[2])
    })
  structure(bg, pval = p.val, pvalue = p.val, df = df, tname = "Breusch-Godfrey", class = "ac_test_")
}

#' @export
ac_test_bg.default <- function(x, order, type, fill) {
  stop("not available method.", call. = FALSE)
}


# Autocorrelation test ----------------------------------------------------


#' Autocorrelation tests
#'
#' @param x the residuals or an `ivx` object.
#' @param lag_max the maximum length of lags.
#'
#' @name ac_test
#' @export
#' @examples
#' obj <- ivx(hpi ~ cpi + def + int + log(res), data = ylpc)
#' if (requireNamespace("lmtest", quietly = TRUE)) {
#'   lmtest::bgtest(hpi ~ cpi + def + int + log(res), data = ylpc)
#' }
#' ac_test(obj, 5)
#'
ac_test <- function(x, lag_max = 5) {
  UseMethod("ac_test")
}

#' @export
ac_test.ivx <- function(x, lag_max = 5) {
  res <- ols_residuals(x)
  stats <- ac_test.default(res, lag_max)
  bg <- pval <-  vector("numeric", lag_max)
  for(i in 1:lag_max) {
    temp <- ac_test_bg.ivx(x, order = i)
    pval[i] <- attr(temp, "pval")
    bg[i] <- unclass(temp)
  }
  stats$BreuschGodfrey <- bg
  attr(stats, "pvalue")$BreuschGodfrey  <- pval
  stats
}

#' @export
ac_test.default <- function(x, lag_max = 5) {
  bp <- lb <- vector("numeric", lag_max)
  for(i in 1:lag_max) {
    lb[i] <- Box.test(x, lag = i, type = "Ljung-Box")$statistic
    bp[i] <- Box.test(x, lag = i, type = "Box-Pierce")$statistic
  }
  stats <- data.frame(
    Lag = 1:lag_max,
    Wald = unclass(ac_test_wald(x, 1:lag_max)),
    LjungBox = lb,
    BoxPierce = bp
  )
  pval <- data.frame(
    Lag = 1:lag_max,
    Wald = 1 - pchisq(stats[,2], 1:lag_max),
    LjungBox = 1 - pchisq(stats[,3], 1:lag_max),
    BoxPierce = 1 - pchisq(stats[,4], 1:lag_max)
  )
  structure(
    stats,
    pvalue = pval,
    class = c("ac_test", "data.frame")
  )
}

stars_pval <- function (pval) {
  unclass(
    symnum(pval, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.01, 0.05, 0.1, 1),
           symbols = c("***", "**", "*", " "))
  )
}

#' @export
print.ac_test <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  pval <- attr(x, "pvalue")
  out <- data.frame(Lag = x$Lag)
  for (nm in setdiff(names(x), "Lag")) {
    out[[nm]] <- paste0(formatC(x[[nm]], digits = digits), stars_pval(pval[[nm]]))
  }
  print(out, row.names = FALSE)
}

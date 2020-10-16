
# AR Wald Test -------------------------------------------------------------


#' Autocorrelation test from yplc
#'
#'
#'
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
  war <- war[lag]
  names(war) <- lag
  attr(war, "pval") <- 1 - pchisq(war, lag)
  class(war) <- "ac_test_wald"
  war
}

#' @export
ac_test_wald.ivx <- function(x, lag = 1) {
  res_ols <- x$residuals_ols
  ac_test_wald.default(x, lag = lag)
}


# Ljung-Box ---------------------------------------------------------------

ac_test_lb <- function(x, lag) {
  UseMethod("ac_test_lb")
}

ac_test_lb.default <- function(x, lag) {
  max_lag <- max(lag)
  lb <- vector("numeric", max_lag)
  for(i in 1:max_lag) {
    lb[i] <- Box.test(x, lag = i, type = "Ljung-Box")$statistic
  }
  out <- lb[lag]
  names(out) <- lag
  attr(out, "pval") <-  1 - pchisq(lb[,3], 1:max)
  class(out) <- "ac_test_lb"
  out
}

ac_test_lb.ivx <- function(x, lag) {
  res <- x$residuals_ols
  ac_test_lb.default(res, lag)
}

print.ac_test_lb <- function(x, digits = max(3L, getOption("digits") - 3L)) {
  pval <- attr(x, "pvalue")
  out <- paste0(formatC(x[,i], digits = digits), stars)
  print(out, row.names = FALSE)
}

# Box-Pierce --------------------------------------------------------------


ac_test_bp <- function(x, lag) {
  UseMethod("ac_test_bp")
}

ac_test_bp.default <- function(x, lag) {
  bp <- vector("numeric", max)
  for(i in 1:max(lag)) {
    bp[i] <- Box.test(res, lag = i, type = "Box-Pierce")$statistic
  }
  out <- bp[lag]
  names(out) <- lag
  attr(out, "pval") <- 1 - pchisq(bp, lag)
  out
}

ac_test_bp.ivx <- function(x, lag) {
  res <- x$residuals_ols
  ac_test_bp.default(x, lag)
}

print.ac_test_bp <- function(x, digits = max(3L, getOption("digits") - 3L)) {
  pval <- attr(x, "pvalue")
  out <- paste0(formatC(x[,i], digits = digits), stars)
  print(out, row.names = FALSE)
}

# Breusch–Godfrey ----------------------------------------------------------

ac_test_bg.ivx <- function(x) {
  X <- model.frame(x)
  Y <- model.matrix(x)
  res <- x$residuals_ols
}


bgtest <- function (x, order = 1, fill = 0) {
  # res <- x$residuals_ols
  y <- x[-1]
  n <- NROW(res)
  Z <- sapply(order, function(x) c(rep(fill, length.out = x), res[1:(n - x)]))
  auxfit <- lm.fit(cbind(X, Z), res)

  bg <- n * sum(auxfit$fitted.values^2)/sum(res^2)
  return(bg)
  # p.val <- pchisq(bg, m, lower.tail = FALSE)
}

# lmtest::bgtest(hpi.ols.resid[2:173] ~ hpi.ols.resid[1:172], order = 1, type = "Chisq")



# Autocorrelation test ----------------------------------------------------

#' If p-value < 0.051: You can reject the null hypothesis assuming a 5% chance of making a mistake.
#' So you can assume that your values are showing dependence on each other.
#'
#' x <- ivx(hpi ~ cpi + def + int + log(res), data = ylpc)
#' lmtest::bgtest(hpi ~ cpi + def + int + log(res), data = ylpc)
#' ac_test(x, 5)
#' ac_test.ivx(x)
#' t2 <- ivx_ar(hpi ~ cpi + def + int + res, data = ylpc)
#'
#' @export
ac_test <- function(x, max = 5) {
  UseMethod("ac_test")
}

#' @export
ac_test.ivx <- function(x, max = 5) {
  res <- x$residuals_ols
  ac_test.default(res, max)
}

#' @export
ac_test.default <- function(x, max = 5) {
  lb <- vector("numeric", max)
  bp <- vector("numeric", max)
  for(i in 1:max) {
    lb[i] <- Box.test(x, lag = i, type = "Ljung-Box")$statistic
    bp[i] <- Box.test(x, lag = i, type = "Box-Pierce")$statistic
  }
  stats <- data.frame(
    Lag = 1:max,
    Wald = ac_test_wald(x, 1:max),
    LjungBox = lb,
    BoxPierce = bp
  )
  pval <- data.frame(
    Lag = 1:max,
    Wald = 1 - pchisq(stats[,2], 1:max),
    LjungBox = 1 - pchisq(stats[,3], 1:max),
    BoxPierce = 1 - pchisq(stats[,4], 1:max)
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
print.ac_test <- function(x, digits = max(3L, getOption("digits") - 3L)) {
  pval <- attr(x, "pvalue")
  lst <- list()
  for(i in 2:4) {
    stars <- stars_pval(pval[,i])
    lst[[i]] <- paste0(formatC(x[,i], digits = digits), stars)
  }
  out <- data.frame(
    Lag = 1:nrow(x),
    Wald = lst[[2]],
    LjungBox = lst[[3]],
    BoxPierce = lst[[4]]
  )
  print(out, row.names = FALSE)
}

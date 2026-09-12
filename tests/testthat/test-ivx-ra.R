context("test-ivx-ra")

test_that("ivx_ra fits and the ivx methods work on it", {
  m <- ivx_ra(Ret ~ DP, data = kms)
  expect_s3_class(m, c("ivx_ra", "ivx"))
  expect_equal(m$Wald_Joint, unname(m$tstat^2))
  expect_true(m$ar_order >= 1 && m$ar_order <= 5)
  expect_error(capture.output(print(m)), NA)
  s <- summary(m)
  expect_equal(colnames(coef(s)), c("Estimate", "Std. Error", "t value", "Wald Ind", "Pr(> chi)"))
  expect_equal(dim(vcov(m)), c(1, 1))
  expect_length(delta(m), 1)

  m2 <- ivx_ra(Ret ~ DP + TBL, data = kms, ar = 2, ar_ic = "bic")
  expect_equal(m2$ar_order, 2)
  expect_equal(m2$ar_method, "fixed")
  expect_length(coef(m2), 2)
  expect_length(m2$gamma, 2)
  expect_equal(dim(vcov(m2)), c(2, 2))
  expect_true(m2$Wald_Joint > 0)
})

test_that("ivx_ra validates ar and matches ivx_ra_fit", {
  expect_error(ivx_ra(Ret ~ DP, data = kms, ar = 0), "positive integer")
  expect_error(ivx_ra(Ret ~ DP, data = kms, ar = 1.5), "positive integer")
  f <- ivx_ra_fit(kms$Ret, as.matrix(kms[, c("DP", "TBL")]), ar = 3)
  m <- ivx_ra(Ret ~ DP + TBL, data = kms, ar = 3)
  expect_equal(unname(coef(m)), unname(f$coefficients))
})

test_that("residual augmentation removes the endogeneity bias (true innovations)", {
  # with the true innovations partialled out the estimator is unbiased; check the
  # mechanics on one long simulated sample: gamma ~ -0.95 and |T beta| small
  set.seed(11)
  Tn <- 5000
  e <- rnorm(Tn); u <- -0.95 * e + sqrt(1 - 0.95^2) * rnorm(Tn)
  x <- cumsum(e); y <- c(0, u[-1])
  f <- ivx_ra_fit(y, matrix(x), ar = 1)
  expect_equal(unname(f$gamma), -0.95, tolerance = 0.02)
  expect_lt(abs(f$tstat), 3)
})

test_that("horizon > 1 is the transformed-regression estimator of DRT (2023)", {
  # trf_sum is A_h' applied to z restricted to t <= n - h + 1
  z <- rnorm(20); h <- 4; n <- 20
  A <- matrix(0, n - h + 1, n)
  for (i in 1:(n - h + 1)) A[i, i:(i + h - 1)] <- 1
  expect_equal(trf_sum(z, h), drop(t(A) %*% z[1:(n - h + 1)]))
  expect_equal(trf_sum(z, 1), z)

  # scalar eq. (4.9) written out directly
  x <- kms$DP; y <- kms$Ret; n <- length(y); h <- 6
  rho <- 1 - 1 / (n - 1)^0.95
  zf <- c(0, stats::filter(diff(x), rho, "recursive"))
  eps <- lm.fit(cbind(x[-n]), x[-1])$residuals
  eps <- eps - mean(eps)
  idx <- 2:n; nn <- length(idx)
  yd <- y[idx] - mean(y[idx]); xd <- x[idx - 1] - mean(x[idx - 1]); zl <- zf[idx - 1]
  ytil <- yd - eps * lm.fit(cbind(eps), yd)$coefficients
  ztr <- sapply(seq_len(nn), function(t) sum(zl[max(1, t - h + 1):min(t, nn - h + 1)]))
  b <- sum(ztr * ytil) / sum((zl * xd)[1:(nn - h + 1)])
  f <- ivx_ra_fit(y, matrix(x), ar = 1, horizon = h)
  expect_equal(unname(coef(f)), b)
  expect_equal(f$horizon, h)

  # h = 1 is unchanged
  expect_equal(coef(ivx_ra(Ret ~ DP, data = kms, ar = 1, horizon = 1)),
               coef(ivx_ra(Ret ~ DP, data = kms, ar = 1)))
  expect_error(ivx_ra(Ret ~ DP, data = kms, horizon = 0), "positive integer")
})

context("test-elliott-cf")

test_that("elliott_cf fits the augmented regression and reports the Wald test", {
  m <- elliott_cf(Ret ~ DP + TBL, ~ LTY, data = kms, lags = 1)
  expect_s3_class(m, "elliott_cf")
  expect_named(m$coefficients, c("DP", "TBL"))
  expect_named(m$gamma, c("LTY", "LTY_lag1"))
  expect_equal(m$df, 2)
  expect_equal(m$p.value, 1 - pchisq(m$Wald, 2))
  expect_true(all(abs(m$rho_resid) <= 1))
  expect_error(capture.output(print(m)), NA)
  # coefficients are OLS of y_t on x_{t-1}, z_t, z_{t-1}
  n <- nrow(kms); idx <- 3:n
  f <- lm(kms$Ret[idx] ~ kms$DP[idx - 1] + kms$TBL[idx - 1] + kms$LTY[idx] + kms$LTY[idx - 1])
  expect_equal(unname(m$coefficients), unname(coef(f)[2:3]))
  expect_equal(unname(m$gamma), unname(coef(f)[4:5]))
  # non-robust standard errors are the lm ones
  m0 <- elliott_cf(Ret ~ DP, ~ TBL, data = kms, robust = FALSE)
  f0 <- lm(kms$Ret[-1] ~ kms$DP[-n] + kms$TBL[-1])
  expect_equal(unname(m0$se), unname(summary(f0)$coefficients[2, 2]))
  expect_error(elliott_cf_fit(kms$Ret, kms$DP, kms$TBL, lags = -1), "lags")
})

test_that("an orthogonalising covariate restores the chi-square null", {
  set.seed(8); n <- 300
  v <- rnorm(n); u <- -0.9 * v + sqrt(0.19) * rnorm(n)
  x <- cumsum(v); y <- u
  m <- elliott_cf_fit(y, x, v)          # z_t = v_t orthogonalises exactly
  expect_lt(abs(m$rho_resid), 0.15)
  expect_equal(unname(m$gamma), -0.9, tolerance = 0.1)
  # a noisy proxy leaves correlation behind, which the diagnostic reports
  expect_gt(abs(elliott_cf_fit(y, x, v + 0.3 * rnorm(n))$rho_resid), 0.3)
})

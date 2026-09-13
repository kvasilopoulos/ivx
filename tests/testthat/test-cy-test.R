context("test-cy-test")

test_that("cy_test runs and its pieces are consistent", {
  m <- cy_test(Ret ~ DP, data = kms)
  expect_s3_class(m, "cy_test")
  expect_true(m$ci["lower"] < m$ci["upper"])
  expect_true(m$c_ci["lower"] < m$c_ci["upper"])
  expect_equal(unname(m$rho_ci), 1 + unname(m$c_ci) / nrow(kms))
  expect_equal(unname(m$reject), unname(c(m$ci["lower"] > 0, m$ci["upper"] < 0)))
  expect_true(abs(m$delta) <= 1)
  expect_error(capture.output(print(m)), NA)
  expect_error(cy_test(Ret ~ DP + TBL, data = kms), "single predictor")
  # sign flip: negating the predictor negates beta and the interval
  f <- cy_test_fit(kms$Ret, -kms$DP)
  expect_equal(unname(f$ci), -rev(unname(m$ci)), tolerance = 1e-8)
  expect_true(f$flipped)
})

test_that("DF-GLS statistic and its inversion behave", {
  # 5% asymptotic critical value at c = 0 is about -1.95 (ERS, 1996)
  expect_equal(unname(dfgls_q["0.05", "0"]), -1.95, tolerance = 0.03)
  # inversion is monotone: a larger statistic gives a larger c
  expect_lt(dfgls_invert(-3, 0.5), dfgls_invert(-1, 0.5))
  expect_equal(dfgls_invert(-100, 0.5), -100)
  # Table 2 interpolation
  expect_equal(unname(cy_alpha1(-0.9)), c(0.060, 0.130))
  expect_equal(unname(cy_alpha1(-0.9125)), c(0.0575, 0.1225))
  # ar_var: AR(1) variance s2 / (1 - phi^2)
  expect_equal(ar_var(0.5, 2), 2 / 0.75)
  # DF-GLS on a random walk is O(1), on white noise strongly negative
  set.seed(2); n <- 400
  expect_gt(dfgls_stat(cumsum(rnorm(n)), 1), -4)
  expect_lt(dfgls_stat(rnorm(n), 1), -10)
})

test_that("Q-estimate reduces to OLS when delta = 0 and the rho interval is a point", {
  # with sigma_ue = 0 the correction term vanishes, so beta(rho) is the OLS slope
  set.seed(3); n <- 300
  x <- cumsum(rnorm(n)); y <- rnorm(n)
  f <- cy_test_fit(y, x, lag_max = 1)
  ols <- unname(coef(lm(y[-(1:2)] ~ x[-c(1, n)]))[2])
  # beta(rho) = ols - k * cov(x_t - rho x_{t-1}, x^mu); |k| is small when delta ~ 0
  expect_equal(unname(f$beta_rho[1]), ols, tolerance = 0.2)
  expect_equal(f$estimate, ols)
})

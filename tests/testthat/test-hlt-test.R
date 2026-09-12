context("test-hlt-test")

test_that("hlt_test runs, switches and prints", {
  h <- hlt_test(Ret ~ DP, data = kms)
  expect_s3_class(h, "hlt_test")
  expect_true(h$test %in% c("T_N", "T_con", "T~_con"))
  expect_true(is.finite(h$statistic) && is.finite(h$cv))
  expect_true(abs(h$rho_xy) <= 1)
  expect_true(h$adf_lag >= 0 && h$adf_lag <= floor(12 * (nrow(kms) / 100)^0.25))
  expect_error(capture.output(print(h)), NA)
  hl <- hlt_test(Ret ~ DP, data = kms, alternative = "less", level = 0.1)
  expect_true(hl$cv < 0)
  expect_error(hlt_test(Ret ~ DP + TBL, data = kms), "single predictor")
  expect_error(hlt_test(Ret ~ DP, data = kms, level = 0.2), "level")
})

test_that("statistics match hand computations", {
  x <- kms$DP; y <- kms$Ret; n <- length(x)
  h <- hlt_test_fit(y, x)
  # T is the usual OLS t-ratio
  expect_equal(h$t, unname(summary(lm(y[-1] ~ x[-n]))$coefficients[2, 3]))
  expect_equal(h$estimate, unname(coef(lm(y[-1] ~ x[-n]))[2]))
  # T~: quasi-GLS mean with c_bar = 7, no intercept in the regression
  rb <- 1 - 7 / n
  mu <- lm.fit(cbind(c(1, rep(1 - rb, n - 1))), c(x[1], x[-1] - rb * x[-n]))$coefficients
  xg <- x[-n] - mu; yd <- y[-1] - mean(y[-1])
  expect_equal(h$t_gls, unname(summary(lm(yd ~ xg - 1))$coefficients[1, 3]))
  # response surface at rho = 0 is the intercept
  expect_equal(hlt_cv(0, 0.05, FALSE), 1.707)
  expect_equal(hlt_cv(0, 0.01, TRUE), 2.377)
})

test_that("weak persistence selects the normal test, strong persistence the conservative ones", {
  set.seed(3); n <- 300
  v <- rnorm(n); u <- -0.9 * v + sqrt(0.19) * rnorm(n)
  expect_equal(hlt_test_fit(u, stats::filter(v, 0.2, "recursive"))$test, "T_N")
  expect_equal(hlt_test_fit(u, cumsum(v))$test, "T~_con")
  expect_equal(hlt_test_fit(u, cumsum(v), alternative = "less")$test, "T_con")
  expect_equal(hlt_test_fit(u, cumsum(v), alternative = "less")$cv, -hlt_cv(0.9, 0.05, FALSE), tolerance = 0.1)
})

context("test-el-test")

test_that("el_ratio is the Owen EL ratio", {
  set.seed(1)
  z <- matrix(rnorm(80), 80)
  # -2 log EL for a mean: solve the dual by hand and compare
  f <- function(l) -sum(log(1 + l * z))
  l <- optimize(f, c(-0.5, 0.5))$minimum
  expect_equal(el_ratio(z), 2 * sum(log(1 + l * z)), tolerance = 1e-5)
  expect_equal(el_ratio(z - mean(z)), 0, tolerance = 1e-8)
  # mean far from 0 gives a large statistic
  expect_gt(el_ratio(z + 2), 50)
})

test_that("el_test runs and returns chi-square p-values", {
  e <- el_test(Ret ~ DP, data = kms)
  expect_s3_class(e, "el_test")
  expect_named(e$stat, c("beta2", "beta1", "joint"))
  expect_true(all(e$stat >= 0))
  expect_equal(unname(e$p.value), 1 - pchisq(unname(e$stat), c(1, 1, 2)))
  expect_equal(e$m, floor(nrow(kms) / 2))
  expect_error(capture.output(print(e)), NA)
  expect_error(el_test(Ret ~ DP + TBL, data = kms), "single predictor")
  # profile statistics are at most the joint one
  expect_lte(e$stat[["beta2"]], e$stat[["joint"]] + 1e-8)
  expect_lte(e$stat[["beta1"]], e$stat[["joint"]] + 1e-8)
})

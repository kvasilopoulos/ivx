context("test-ivx-iv")

test_that("ivx_iv fits with each instrument set and the ivx methods work", {
  for (inst in c("comb", "sin", "frac", "diff")) {
    m <- ivx_iv(Ret ~ DP + TBL, data = kms, instruments = inst)
    expect_s3_class(m, c("ivx_iv", "ivx"))
    expect_length(coef(m), 2)
    expect_equal(NCOL(m$instruments), if (inst == "comb") 4 else 2)
    expect_equal(unname(m$Wald_Ind), unname(m$tstat^2))
    expect_true(is.finite(m$Wald_Joint))
    expect_error(capture.output(print(m)), NA)
    expect_error(capture.output(print(summary(m))), NA)
  }
  expect_equal(unname(coef(ivx_iv(Ret ~ DP, data = kms))),
               unname(ivx_iv_fit(kms$Ret, as.matrix(kms$DP))$coefficients))
  expect_error(ivx_iv(Ret ~ DP, data = kms, d = 0.7), "d")
})

test_that("instrument constructions are right", {
  set.seed(1); x <- matrix(rnorm(40), 20); d <- 0.5; n <- 20
  pi <- cumprod(c(1, (seq_len(n - 1) - 1 - d) / seq_len(n - 1)))
  ref <- sapply(1:2, function(j) sapply(1:n, function(t) sum(pi[1:t] * x[t:1, j])))
  expect_equal(frac_diff(x, d), ref)
  expect_equal(frac_diff(matrix(1:5), 1)[, 1], rep(1, 5))       # d = 1: first difference
  expect_equal(long_diff(matrix(1:10), 0.2, 0.85)[, 1], rep(1, 10))  # k = 1
  expect_equal(long_diff(matrix(1:10), 3, 1)[, 1], 1:10)  # k = 30 > t - 1: x_t - x_0
  expect_equal(dim(sin_inst(30, 2)), c(30, 2))
})

test_that("2SLS with a single instrument equals the IV formula", {
  x <- kms$DP; y <- kms$Ret; n <- length(x)
  m <- ivx_iv_fit(y, matrix(x), instruments = "sin")
  z <- sin(pi * seq_len(n - 1) / (n - 1))
  xd <- x[-n] - mean(x[-n]); yd <- y[-1] - mean(y[-1]); zd <- z - mean(z)
  expect_equal(unname(coef(m)), sum(zd * yd) / sum(zd * xd))
})

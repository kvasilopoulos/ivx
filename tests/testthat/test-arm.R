context("test-arm")

test_that("arm fits and the ivx methods work on it", {
  m <- arm(Ret ~ DP + TBL, data = kms)
  expect_s3_class(m, c("arm", "ivx"))
  expect_length(coef(m), 2)
  expect_length(m$phi, 2)
  expect_equal(dim(m$Phi), c(2, 2))
  expect_true(all(Mod(eigen(m$Phi)$values) < 1))
  expect_equal(unname(m$Wald_Ind), unname(m$tstat^2))
  expect_error(capture.output(print(m)), NA)
  s <- summary(m)
  expect_equal(colnames(coef(s)), c("Estimate", "Std. Error", "t value", "Wald Ind", "Pr(> chi)"))
  expect_equal(dim(vcov(m)), c(2, 2))
  expect_equal(unname(coef(m)), unname(arm_fit(kms$Ret, as.matrix(kms[, c("DP", "TBL")]))$coefficients))
})

test_that("Nicholls-Pope bias reduces to Kendall's 1 + 3 phi and the correction is applied", {
  phi <- 0.9; s2 <- 2
  expect_equal(np_bias(matrix(phi), matrix(s2)), matrix(1 + 3 * phi))
  f <- arm_fit(kms$Ret, as.matrix(kms$DP), iter = 1)
  n <- nrow(kms) - 1
  expect_equal(drop(f$Phi), drop(f$Phi_ols) + (1 + 3 * drop(f$Phi_ols)) / n, tolerance = 1e-6)
  # single predictor: beta_c is the OLS slope of y on (x_{t-1}, v_c)
  x <- kms$DP; y <- kms$Ret; nx <- length(x)
  vc <- x[-1] - (mean(x[-1]) - drop(f$Phi) * mean(x[-nx])) - drop(f$Phi) * x[-nx]
  b <- lm.fit(cbind(1, x[-nx], vc), y[-1])$coefficients
  expect_equal(unname(coef(f)), unname(b[2]))
  expect_equal(unname(f$phi), unname(b[3]))
})

test_that("arm removes most of the Stambaugh bias", {
  set.seed(4)
  n <- 100; R <- 200
  b <- replicate(R, {
    v <- rnorm(n); u <- -0.9 * v + sqrt(1 - 0.81) * rnorm(n)
    x <- stats::filter(v, 0.9, "recursive"); y <- u
    c(arm = unname(arm_fit(y, matrix(x))$coefficients),
      ols = unname(lm.fit(cbind(1, x[-n]), y[-1])$coefficients[2]))
  })
  expect_lt(abs(mean(b["arm", ])), abs(mean(b["ols", ])) / 3)
})

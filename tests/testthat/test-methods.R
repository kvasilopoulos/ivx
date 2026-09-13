context("test-methods")

obj <- ivx(Ret ~ TBL + EP + LTY, kms)
obj_ar <- ivx_ar(Ret ~ DP + LTY, kms)

test_that("summary.ivx return the same with ivx", {
  expect_equal(residuals(summary(obj)), residuals(obj))
  expect_equal(fitted(summary(obj)), fitted(obj))

  expect_equal(residuals(summary(obj_ar)), residuals(obj_ar))
  expect_equal(fitted(summary(obj_ar)), fitted(obj_ar))

})


test_that("step-methods",{
  expect_error(drop1(obj), NA)
  expect_error(add1(obj, "DE"), NA)
  expect_error(capture.output(step(obj)), NA)
  expect_error(deviance(obj), NA)
  expect_error(logLik(obj), NA)
})


test_that("lm-style methods work on ivx fits", {
  expect_equal(length(case.names(obj)), length(residuals(obj)))
  expect_equal(dim(model.matrix(obj)), c(nrow(kms), 3))
  expect_equal(nrow(model.frame(obj)), nrow(kms))
  expect_true(is.finite(extractAIC(obj)[2]))
  expect_equal(unname(deviance(obj)), sum(residuals(obj)^2))
  expect_error(capture.output(print(summary(obj_ar))), NA)
  expect_error(capture.output(print(obj_ar)), NA)
  expect_error(capture.output(print(delta(obj))), NA)
  expect_equal(dim(vcov(summary(obj))), c(3, 3))
})

test_that("weighted fits handle zero weights and offsets", {
  w <- rep(1, nrow(kms)); w[1:5] <- 0
  m <- ivx(Ret ~ DP + TBL, data = kms, weights = w)
  expect_length(residuals(m), nrow(kms))
  expect_equal(sum(is.na(residuals(m))), 1)   # first kept observation is lost to the lag
  expect_equal(length(m$weights), nrow(kms))
  expect_equal(unname(coef(m)), unname(coef(ivx(Ret ~ DP + TBL, data = kms[-(1:5), ]))))
  m2 <- ivx(Ret ~ DP, data = kms, offset = rep(0.01, nrow(kms)))
  expect_equal(unname(coef(m2)), unname(coef(ivx(Ret ~ DP, data = kms))), tolerance = 1e-8)
  expect_warning(ivx(Ret ~ DP - 1, data = kms), "intercept")
  expect_error(ivx(cbind(Ret, DP) ~ TBL, data = kms), "multivariate")
})

test_that("texreg extract method builds a texreg object", {
  skip_if_not_installed("texreg")
  tr <- extract.ivx(obj, include.aic = TRUE, include.bic = TRUE, include.rsquared = TRUE, include.adjrs = TRUE)
  expect_s4_class(tr, "texreg")
  expect_equal(length(tr@coef), 3)
})

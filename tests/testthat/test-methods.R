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


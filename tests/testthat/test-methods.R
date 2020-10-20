context("test-methods")

object <- ivx(Ret ~ TBL + EP + LTY, kms)
obj_ar <- ivx(Ret ~ DP + LTY, kms)

test_that("summary.ivx return the same with ivx", {
  expect_equal(residuals(summary(obj)), residuals(obj))
  expect_equal(fitted(summary(obj)), fitted(obj))

  expect_equal(residuals(summary(obj_ar)), residuals(obj_ar))
  expect_equal(fitted(summary(obj_ar)), fitted(obj_ar))

})


test_that("step-methods",{
  expect_error(drop1(object), NA)
  expect_error(add1(object, "DE"), NA)
  expect_error(capture.output(step(object)), NA)
  expect_error(deviance(object), NA)
  expect_error(logLik(object), NA)
})


context("test-warn-stop")

test_that("interecept warning", {
  expect_warning(ivx(Ret ~ - 1 + LTY, data = kms),
                 "ivx estimation does not include an intercept by construction")
})

test_that("empty-colnames",{
  monthly_mat <- as.matrix(kms[, -1]) %>% unname()
  reg <- ivx_fit(kms[,1],  as.matrix(kms[,c(2,3)]))
  expect_equal(c("DE", "LTY"), reg$cnames)
})


test_that("printing coverage", {
  reg <- ivx(Ret ~  LTY, data = monthly)
  expect_error(capture.output(print(reg)), NA)
  expect_error(capture.output(print(summary(reg))), NA)
  expect_error(capture.output(print(residuals(reg))), NA)
  expect_error(capture.output(print(vcov(reg))), NA)
})

test_that("right delta object", {
  reg_lm <- lm(Ret ~  LTY, data = monthly)
  expect_error(delta(reg_lm), "Wrong object")
})

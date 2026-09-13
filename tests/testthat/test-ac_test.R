context("test-ac_test")

obj <- ivx(hpi ~ cpi + inv, data = ylpc)

test_that("autocorrelation tests run on ivx objects and on residual vectors", {
  w <- ac_test_wald(obj, lag = 1:3)
  expect_s3_class(w, "ac_test_")
  expect_length(w, 3)
  expect_true(all(attr(w, "pval") >= 0 & attr(w, "pval") <= 1))
  expect_equal(unname(w), unname(ac_test_wald(obj$ols$residuals, lag = 1:3)))
  for (f in list(ac_test_lb, ac_test_bp)) {
    s <- f(obj, lag = c(1, 3))
    expect_length(s, 2)
    expect_equal(unname(attr(s, "pval")), 1 - pchisq(as.numeric(s), c(1, 3)))
    expect_error(capture.output(print(s)), NA)
  }
  bg <- ac_test_bg(obj, order = 2)
  expect_true(attr(bg, "pval") >= 0 && attr(bg, "pval") <= 1)
  expect_error(ac_test_bg(obj$ols$residuals, order = 2), "not available")
  a <- ac_test(obj, lag_max = 3)
  expect_error(capture.output(print(a)), NA)
})

test_that("Ljung-Box and Box-Pierce match Box.test", {
  r <- obj$ols$residuals
  expect_equal(unname(unclass(ac_test_lb(r, lag = 2)))[1], unname(Box.test(r, 2, type = "Ljung-Box")$statistic))
  expect_equal(unname(unclass(ac_test_bp(r, lag = 2)))[1], unname(Box.test(r, 2, type = "Box-Pierce")$statistic))
})

test_that("ac_test rejects the extensions that do not store OLS residuals", {
  for (m in list(ivx_ar(hpi ~ cpi + inv, data = ylpc), arm(hpi ~ cpi + inv, data = ylpc))) {
    expect_error(ac_test(m), "fits only")
    expect_error(ac_test_bg(m), "fits only")
  }
})

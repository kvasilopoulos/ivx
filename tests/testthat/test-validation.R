context("test-validation")

test_that("fitter functions validate their inputs", {
  y <- kms$Ret; x <- as.matrix(kms$DP)
  expect_error(ivx_fit(y, x[-1, , drop = FALSE]), "incompatible")
  expect_error(ivx_fit(y[1:0], x[0, , drop = FALSE]), "non-NA")
  expect_equal(ivx_fit(y, x[, 0, drop = FALSE])$df.residuals, length(y))
  expect_error(ivx_wfit(y, x, rep(1, 10)), "incompatible")
  expect_error(ivx_wfit(y, x, rep(-1, length(y))), "negative")
  expect_error(ivx(Ret ~ DP, data = kms, weights = rep("a", nrow(kms))), "numeric")
  expect_error(ivx_ra_fit(y, x[-1, , drop = FALSE]), "incompatible")
  expect_error(ivx_sys_fit(cbind(y), x[-1, , drop = FALSE]), "incompatible")
  expect_error(arm_fit(y, x[-1, , drop = FALSE]), "incompatible")
  expect_error(el_test_fit(y[1:10], x[1:10]), "too few")
  expect_error(el_test_fit(y, x[-1]), "incompatible")
  expect_error(cy_test_fit(y, x[-1]), "incompatible")
  expect_error(hlt_test_fit(y, x[-1]), "incompatible")
  expect_error(elliott_cf_fit(y, x, kms$TBL[-1]), "incompatible")
  expect_error(ivx_iv_fit(y, x[-1, , drop = FALSE]), "incompatible")
  expect_error(ivx_qr_fit(y, x, tau = c(0.1, 0.5)), "single value")
})

test_that("ivx_ar validates its arguments and handles ar = 0", {
  expect_error(ivx_ar(Ret ~ DP, kms, ar = "x"), "non-negative integer")
  expect_error(ivx_ar(Ret ~ DP, kms, ar_max = 0), "positive integer")
  expect_error(ivx_ar(Ret ~ DP, kms, ar_grid = 1), "function")
  expect_message(m0 <- ivx_ar(Ret ~ DP, kms, ar = 0), "Using `ivx`")
  expect_s3_class(m0, "ivx")
  expect_equal(coef(m0), coef(ivx(Ret ~ DP, kms)))
  m <- ivx_ar(Ret ~ DP, kms, ar = 2, robust = TRUE)
  expect_error(capture.output(print(summary(m))), NA)
})

test_that("empty model and Breusch-Godfrey F variant", {
  m <- ivx(Ret ~ 1, data = kms)
  expect_length(coef(m), 0)
  expect_error(capture.output(print(m)), NA)
  obj <- ivx(hpi ~ cpi + inv, data = ylpc)
  f <- ac_test_bg(obj, order = 2, type = "F")
  expect_equal(length(attr(f, "df")), 2)
  expect_true(attr(f, "pval") >= 0 && attr(f, "pval") <= 1)
})

test_that("arm falls back to Yule-Walker when the OLS VAR is explosive", {
  set.seed(1); n <- 60
  x <- cumsum(cumsum(rnorm(n)))      # near-I(2): OLS AR(1) root can exceed one
  y <- rnorm(n)
  expect_error(f <- arm_fit(y, matrix(x)), NA)
  expect_true(is.finite(f$tstat))
})

test_that("logLik and extractAIC handle weights", {
  w <- rep(1, nrow(kms)); w[1:3] <- 0
  m <- ivx(Ret ~ DP, data = kms, weights = w)
  ll <- logLik(m)
  expect_true(is.finite(ll))
  expect_equal(attr(ll, "nobs"), nrow(kms) - 3 - 1)
  expect_equal(unname(logLik(ivx(Ret ~ DP, data = kms))), unname(logLik(ivx(Ret ~ DP, data = kms, weights = rep(1, nrow(kms))))))
})


# βOLS tOLSˆ βIVX-KMS Wβ,IVX-KMSˆ βIVX-AR Wβ,IVX-AR θ q Wφ
# CPI −0.0001 −5.211∗∗∗ −0.0001 36.711∗∗∗ −0.0001 1.803 −0.1685 5 134.042∗∗∗
# INV 0.0070 8.750∗∗∗ 0.0078 92.486∗∗∗ 0.0070 13.785∗∗∗ 0.3888 1 17.556∗∗∗

# YLPC Table 1 p6
test_that("ivx_ar univariate", {
  spec <- ivx(hpi ~ inv, ylpc)
  spec_ar <- ivx_ar(hpi ~ inv, ylpc, ar = "forecast")

  expect_equal(unname(spec$coefficients_ols[2]), 0.0070, tol = 0.0001)
  expect_equal(spec$tstat_ols[2,1], 8.750, tol = 0.01)
  expect_equal(unname(spec$coefficients), 0.0078, tol = 0.0001)
  expect_equal(unname(spec$Wald_Ind), 92.486, tol = 0.001)
  expect_equal(delta(spec),  0.388, tol = 0.001)
  expect_equal(unname(spec_ar$coefficients), 0.0070,tol = 0.001)
  expect_equal(unname(spec_ar$Wald_Ind), 13.785,tol = 0.001)
  expect_equal(unname(spec_ar$Wald_AR[1]), 17.566,tol = 0.001)

})

test_that("ar estimation is identical", {
  ar_fixed <- ivx_ar(hpi ~ inv, ylpc, ar = 4)$coefficients_ar
  ar_forecast <- ivx_ar(hpi ~ inv, ylpc, ar = "forecast",
                       stepwise=FALSE, approximation  = FALSE)$coefficients_ar
  ar_auto <- ivx_ar(hpi ~ inv, ylpc, ar = "auto")$coefficients_ar
  expect_equal(ar_fixed, ar_auto)
  expect_equal(ar_forecast, ar_auto)
})

# test_that("ar estimation routines give the same result", {
#   mdl <- ivx(hpi ~ inv, ylpc)
#   res <- mdl$residuals_ols
#   arima(res, order = c(4,0,0), include.mean = FALSE)
#   auto_ar(res)
#   forecast::auto.arima(res, max.q = 0, d = 0, ic = "bic", stepwise=FALSE, approximation  = FALSE)
# })


# test_that("ivx_ar multivariate", {
#   spec <- ivx_ar(hpi ~ inv + une, ylpc, ar =1)
#
#   spec_ar <- ivx_ar(hpi ~ cpi, ylpc, ar = 5)
#   spec_ar$coefficients
#   spec_ar$Wald_Ind
#   spec_ar$Wald_AR
#
# })


# YLPC table6 p17
test_that("ac-test", {
  ac <- ac_test(ivx(hpi ~ cpi + def + int + log(res), data =ylpc))
  expect_equal(ac$Wald[1:5], c(27.768, 48.649, 47.696, 126.944, 125.417), tol = 0.001)
  expect_equal(ac$LjungBox[1:5], c(61.855, 91.561, 145.191, 206.823, 236.423), tol = 0.001)
  expect_equal(ac$BoxPierce[1:5], c(60.795, 89.822, 141.919, 201.439, 229.854), tol = 0.001)
})



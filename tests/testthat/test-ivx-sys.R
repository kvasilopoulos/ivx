context("test-ivx-sys")

test_that("ivx_sys reduces to ivx for a single response, at short and long horizon", {
  for (h in c(1, 4)) {
    u <- ivx_sys(Ret ~ DP + TBL, data = kms, horizon = h)
    m <- ivx(Ret ~ DP + TBL, data = kms, horizon = h)
    expect_equal(drop(coef(u)), coef(m))
    expect_equal(u$Wald_Joint, m$Wald_Joint)
    expect_equal(unname(u$Wald_Eq), m$Wald_Joint)
    expect_equal(unname(vcov(u)), unname(vcov(m)))
    expect_equal(unname(drop(u$delta)), unname(delta(m)))
  }
})

test_that("ivx_sys handles a two-equation system", {
  s <- ivx_sys(cbind(Ret, DE) ~ DP + TBL, data = kms)
  expect_s3_class(s, "ivx_sys")
  expect_equal(dim(coef(s)), c(2, 2))
  expect_equal(dimnames(coef(s)), list(c("Ret", "DE"), c("DP", "TBL")))
  expect_equal(rownames(vcov(s)), c("Ret:DP", "DE:DP", "Ret:TBL", "DE:TBL"))
  expect_equal(s$df, 4)
  expect_named(s$Wald_Eq, c("Ret", "DE"))
  # first equation's Wald equals the univariate joint Wald
  expect_equal(unname(s$Wald_Eq["Ret"]), ivx(Ret ~ DP + TBL, kms)$Wald_Joint)
  expect_error(capture.output(print(s)), NA)
  sm <- summary(s)
  expect_equal(dim(coef(sm)), c(4, 5))
  expect_error(capture.output(print(sm)), NA)
  f <- ivx_sys_fit(as.matrix(kms[, c("Ret", "DE")]), as.matrix(kms[, c("DP", "TBL")]))
  expect_equal(unname(f$coefficients), unname(coef(s)))
})

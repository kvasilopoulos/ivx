context("test-ivx-qr")

skip_if_not_installed("quantreg")

test_that("ivx_qr fits at one and several quantiles", {
  m <- ivx_qr(Ret ~ DP, data = kms, tau = 0.5)
  expect_s3_class(m, c("ivx_qr", "ivx"))
  expect_equal(m$tau, 0.5)
  expect_equal(m$Wald_Joint, unname(m$tstat^2))
  expect_true(m$sparsity > 0)
  expect_true(abs(m$rho_tau) <= 1)
  expect_s3_class(m$rq, "rq")
  expect_error(capture.output(print(m)), NA)
  s <- summary(m)
  expect_s3_class(s, "summary.ivx_qr")
  expect_error(capture.output(print(s)), NA)
  expect_equal(dim(vcov(m)), c(1, 1))

  ms <- ivx_qr(Ret ~ DP + TBL, data = kms, tau = c(0.1, 0.9))
  expect_named(ms, c("0.1", "0.9"))
  expect_length(coef(ms[["0.9"]]), 2)
  expect_equal(ms[["0.9"]]$call$tau, 0.9)
})

test_that("ivx_qr validates tau and passes tuning", {
  expect_error(ivx_qr(Ret ~ DP, data = kms, tau = 1), "tau")
  m <- ivx_qr(Ret ~ DP, data = kms, tau = 0.5, beta = 0.8, cz = 5)
  expect_equal(m$AR$Rz, 1 - 5 / (nrow(kms) - 1)^0.8)
  f <- ivx_qr_fit(kms$Ret, as.matrix(kms$DP), tau = 0.5, beta = 0.8, cz = 5)
  expect_equal(unname(coef(m)), unname(f$coefficients))
})

test_that("ivx_qr_boot resamples blocks of (y, z) and returns percentile inference", {
  m <- ivx_qr(Ret ~ DP + TBL, data = kms, tau = 0.5)
  bt <- ivx_qr_boot(m, B = 49, seed = 1)
  expect_s3_class(bt, "ivx_qr_boot")
  expect_equal(dim(bt$boot), c(49, 2))
  expect_equal(bt$block, ceiling((nrow(kms) - 1)^(1 / 4)))
  expect_equal(dim(bt$ci), c(2, 2))
  expect_true(all(bt$ci[, 1] <= bt$ci[, 2]))
  expect_true(all(bt$p.value >= 0 & bt$p.value <= 1))
  expect_equal(bt$coefficients, coef(m))
  expect_equal(ivx_qr_boot(m, B = 49, seed = 1)$boot, bt$boot)
  expect_error(capture.output(print(bt)), NA)
  # block = n reproduces the original fit in every replication
  n1 <- nrow(kms) - 1
  b1 <- ivx_qr_boot(m, B = 3, block = n1, seed = 2)
  expect_equal(unname(b1$boot[1, ]), unname(coef(m)), tolerance = 1e-6)
  expect_error(ivx_qr_boot(m, block = 0), "block")
  expect_error(ivx_qr_boot(ivx(Ret ~ DP, data = kms)), "ivx_qr")
})

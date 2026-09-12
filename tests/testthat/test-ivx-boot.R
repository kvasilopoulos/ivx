context("test-ivx-boot")

mod <- ivx(Ret ~ DP + TBL, data = kms)

test_that("rwb and frwb return valid p-values", {
  for (type in c("rwb", "frwb")) {
    b <- ivx_boot(mod, B = 49, type = type, seed = 1)
    expect_s3_class(b, "ivx_boot")
    expect_equal(dim(b$boot$tstat), c(49, 2))
    p <- c(b$p.value$Wald_Joint, b$p.value$Wald_Ind, b$p.value$tstat)
    expect_true(all(p >= 0 & p <= 1))
    expect_equal(b$p.value$tstat[, "less"] + b$p.value$tstat[, "greater"],
                 c(DP = 1, TBL = 1), tolerance = 1 / 49)
    expect_error(capture.output(print(b)), NA)
  }
})

test_that("bootstrap is reproducible with seed and works at long horizon", {
  b1 <- ivx_boot(mod, B = 19, seed = 42)
  b2 <- ivx_boot(mod, B = 19, seed = 42)
  expect_equal(b1$boot, b2$boot)
  expect_error(ivx_boot(ivx(Ret ~ DP, kms, horizon = 12), B = 9, seed = 1), NA)
})

test_that("bootstrap rejects unsupported objects", {
  expect_error(ivx_boot(ivx_ar(Ret ~ DP, kms, ar = 1), B = 9), "plain")
  expect_error(ivx_boot(ivx(Ret ~ DP, kms, weights = rep(1, nrow(kms))), B = 9), "weighted")
})

test_that("parallel bootstrap matches serial layout and is reproducible", {
  skip_on_cran()
  # PSOCK workers need an installed ivx; skip when running from a dev checkout only
  b1 <- tryCatch(ivx_boot(mod, B = 20, seed = 7, cores = 2),
                 error = function(e) skip("parallel workers could not load ivx"))
  b2 <- ivx_boot(mod, B = 20, seed = 7, cores = 2)
  expect_equal(dim(b1$boot$tstat), c(20, 2))
  expect_equal(b1$boot, b2$boot)
})

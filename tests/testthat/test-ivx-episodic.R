context("test-ivx-episodic")

mod <- ivx(Ret ~ DP, data = kms)

test_that("subsample statistic reproduces the full-sample IVX t without correction", {
  x <- model.matrix(mod); y <- model.response(model.frame(mod)); n <- nrow(x)
  rho <- 1 - 1 / (n - 1)^0.95
  z <- ivx:::ivx_instrument(x, rho)
  s <- ivx:::sub_ivx_cpp(y[-1], x[-n, , drop = FALSE], z, 0L, as.integer(n - 2), FALSE)
  yd <- y[-1] - mean(y[-1]); xd <- x[-n] - mean(x[-n])
  b <- sum(z * yd) / sum(z * xd); u <- resid(lm(yd ~ xd))
  expect_equal(s[1, 2], b / (sqrt(mean(u^2) * sum(z^2)) / sum(z * xd)))
  expect_equal(s[1, 1], s[1, 2]^2)
})

test_that("ivx_episodic runs for each scheme and returns valid p-values", {
  for (sc in c("rolling", "forward", "backward")) {
    e <- ivx_episodic(mod, scheme = sc, window = 0.3, B = 19, seed = 1)
    expect_s3_class(e, "ivx_episodic")
    expect_named(e$statistic, c("sup", "inf", "sup_sq"))
    expect_true(all(e$p.value >= 0 & e$p.value <= 1))
    expect_equal(dim(e$boot), c(19, 3))
    expect_true(all(e$sequence$end > e$sequence$start))
    expect_error(capture.output(print(e)), NA)
  }
  e <- ivx_episodic(mod, scheme = "forward", window = 0.3, B = 9, seed = 1)
  expect_equal(unique(e$sequence$start), 2L)
  expect_equal(max(e$sequence$end), nrow(kms))
})

test_that("multiple predictors give a sup-Wald test; rwb and robust work", {
  e <- ivx_episodic(ivx(Ret ~ DP + TBL, kms), window = 0.4, B = 9, seed = 1, type = "rwb", robust = TRUE)
  expect_named(e$statistic, "sup_wald")
  expect_equal(names(e$sequence), c("start", "end", "Wald", "t_DP", "t_TBL"))
  expect_error(ivx_episodic(ivx(Ret ~ DP, kms, horizon = 4), B = 9), "horizon")
  expect_error(ivx_episodic(mod, window = 1), "window")
})

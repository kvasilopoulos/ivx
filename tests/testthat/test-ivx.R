context("test-ivx")

# KMS table6 p1531 (monthly data)
test_that("univariate regression ", {
  ivx_model <- ivx(Ret ~ DE, data = monthly)
  expect_equal(
    coef(ivx_model) %>%
      unname() %>%
      round(4),
    -0.0033
  )
  expect_equal(
    summary(ivx_model) %>%
      coef() %>%
      `[`(1, 4) %>%
      round(3),
    0.393
  )

  expect_equal(
    delta(ivx_model) %>%
      round(3),
    -0.067
  )
})

# KMS table8 p1537 (monthly data)
test_that("multivariate regression", {
  ivx_model <- ivx(Ret ~ DP + TBL, data = monthly) %>%
    summary()

  expect_equal(
    coef(ivx_model) %>%
      `[`(, 1) %>%
      unname() %>%
      round(4),
    c(0.0061, -0.0807)
  )
  expect_equal(
    drop(ivx_model$Wald_Joint) %>%
      round(3),
    3.644
  )
})


# KMS table11 p1544 (monthly data)
test_that("univariate long-horizon regression", {
  col1 <- c(0.138, .005, 0.472, 0.803, 0.422, 0.637)
  ivx_model <- ivx(Ret ~ DE, data = monthly, horizon = 4)
  expect_equal(
    ivx_model %>% summary() %>% coef() %>% `[`(4) %>% round(3),
    col1[1]
  )
  expect_equal(
    ivx_model %>%
      update(horizon = 12) %>%
      summary() %>%
      coef() %>%
      `[`(4) %>%
      round(3),
    col1[2]
  )
  expect_equal(
    ivx_model %>% update(horizon = 24) %>% summary() %>%
      coef() %>% `[`(4) %>% round(3),
    col1[3]
  )
  expect_equal(
    ivx_model %>% update(horizon = 36) %>% summary() %>%
      coef() %>% `[`(4) %>% round(3),
    col1[4]
  )
  expect_equal(
    ivx_model %>% update(horizon = 48) %>% summary() %>%
      coef() %>% `[`(4) %>% round(3),
    col1[5]
  )
  expect_equal(
    ivx_model %>% update(horizon = 60) %>% summary() %>%
      coef() %>% `[`(4) %>% round(3),
    col1[6]
  )
})


# KMS table13 p1547 (monthly data)
test_that("multivariate long-horizon regression", {
  col12 <- data.frame(
    col1 = c(5.778, 6.383, 4.990, 4.599, 4.983, 4.321),
    col2 = c(3.894, 3.166, 2.124, 1.915, 1.441, 1.039),
    wald = c(7.638, 7.614, 5.794, 5.383, 5.660, 4.822)
  )
  ivx_model <- ivx(Ret ~ EP + TBL, data = monthly, horizon = 4)
  expect_equal(
    ivx_model %>% summary() %>% coef() %>% `[`(, 4) %>% unname() %>% round(3),
    col12[1,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 12) %>% summary() %>%
      coef() %>% `[`(,4) %>% unname() %>% round(3),
    col12[2,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 24) %>% summary() %>%
      coef() %>% `[`(,4) %>% unname() %>% round(3),
    col12[3,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 36) %>% summary() %>%
      coef() %>% `[`(,4) %>% unname() %>% round(3),
    col12[4,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 48) %>% summary() %>%
      coef() %>% `[`(,4) %>% unname() %>% round(3),
    col12[5,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 60) %>% summary() %>%
      coef() %>% `[`(,4) %>% unname() %>% round(3),
    col12[6,1:2] %>% as.double()
  )
})


test_that("tuning parameters change the instrument and robust changes the vcov", {
  m <- ivx(Ret ~ DP + TBL, data = kms)
  expect_equal(m$tuning, list(beta = 0.95, cz = 1, bandwidth = floor((nrow(kms) - 1)^(1/3))))
  expect_equal(m$AR$Rz, rep(1 - 1 / (nrow(kms) - 1)^0.95, 2))
  expect_equal(m$tstat^2, m$Wald_Ind)
  expect_equal(colnames(coef(summary(m))),
               c("Estimate", "Std. Error", "t value", "Wald Ind", "Pr(> chi)"))

  m2 <- ivx(Ret ~ DP + TBL, data = kms, beta = 0.9, cz = 5, bandwidth = 10)
  expect_equal(m2$AR$Rz, rep(1 - 5 / (nrow(kms) - 1)^0.9, 2))
  expect_equal(m2$tuning$bandwidth, 10)
  expect_false(isTRUE(all.equal(coef(m), coef(m2))))

  mr <- ivx(Ret ~ DP + TBL, data = kms, robust = TRUE)
  expect_true(mr$robust)
  expect_equal(coef(mr), coef(m))
  expect_false(isTRUE(all.equal(vcov(mr), vcov(m))))
  expect_error(ivx(Ret ~ DP, kms, horizon = 4, robust = TRUE), "horizon = 1")
})

test_that("weighted fit honours horizon", {
  w <- ivx(Ret ~ DP + TBL, kms, weights = rep(1, nrow(kms)), horizon = 4)
  expect_equal(coef(w), coef(ivx(Ret ~ DP + TBL, kms, horizon = 4)))
})

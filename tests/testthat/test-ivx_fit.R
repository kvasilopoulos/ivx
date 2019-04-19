context("test-ivx_fit")

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
      `[`(1, 2) %>%
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
    ivx_model %>% summary() %>% coef() %>% `[`(2) %>% round(3),
    col1[1]
  )
  expect_equal(
    ivx_model %>%
      update(horizon = 12) %>%
      summary() %>%
      coef() %>%
      `[`(2) %>%
      round(3),
    col1[2]
  )
  expect_equal(
    ivx_model %>% update(horizon = 24) %>% summary() %>%
      coef() %>% `[`(2) %>% round(3),
    col1[3]
  )
  expect_equal(
    ivx_model %>% update(horizon = 36) %>% summary() %>%
      coef() %>% `[`(2) %>% round(3),
    col1[4]
  )
  expect_equal(
    ivx_model %>% update(horizon = 48) %>% summary() %>%
      coef() %>% `[`(2) %>% round(3),
    col1[5]
  )
  expect_equal(
    ivx_model %>% update(horizon = 60) %>% summary() %>%
      coef() %>% `[`(2) %>% round(3),
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
    ivx_model %>% summary() %>% coef() %>% `[`(, 2) %>% unname() %>% round(3),
    col12[1,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 12) %>% summary() %>%
      coef() %>% `[`(,2) %>% unname() %>% round(3),
    col12[2,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 24) %>% summary() %>%
      coef() %>% `[`(,2) %>% unname() %>% round(3),
    col12[3,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 36) %>% summary() %>%
      coef() %>% `[`(,2) %>% unname() %>% round(3),
    col12[4,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 48) %>% summary() %>%
      coef() %>% `[`(,2) %>% unname() %>% round(3),
    col12[5,1:2] %>% as.double()
  )
  expect_equal(
    ivx_model %>% update(horizon = 60) %>% summary() %>%
      coef() %>% `[`(,2) %>% unname() %>% round(3),
    col12[6,1:2] %>% as.double()
  )
})

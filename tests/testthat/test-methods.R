context("test-methods")

test_that("summary.ivx return the same with ivx", {
  expect_equal(
    ivx(Ret ~ DP + LTY, monthly) %>% summary() %>% residuals(),
    ivx(Ret ~ DP + LTY, monthly) %>% residuals()
  )

  expect_equal(
    ivx(Ret ~ DP + LTY, monthly) %>% summary() %>% fitted(),
    ivx(Ret ~ DP + LTY, monthly) %>% fitted()
  )
})

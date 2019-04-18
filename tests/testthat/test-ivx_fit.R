context("test-ivx_fit")

test_that("Dividend payout ratio", {
  expect_equal(
    ivx(Ret ~ DE, data = monthly) %>% coef() %>% unname() %>% round(4),
    -0.0033
  )
  expect_equal(
    ivx(Ret ~ DE, data = monthly) %>% delta() %>% round(3),
    -0.067
  )

})

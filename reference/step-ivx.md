# Deprecated single-term selection for IVX fits

[`drop1()`](https://rdrr.io/r/stats/add1.html),
[`add1()`](https://rdrr.io/r/stats/add1.html) and hence
[`step()`](https://rdrr.io/r/stats/step.html) on an `ivx` object compare
residual sums of squares of refitted IVX models. IVX is an
instrumental-variable estimator: its residual sum of squares is not
minimised by the fitted coefficients, so RSS-based AIC/F/chi-square
comparisons have no justification and the methods are deprecated. Use
the IVX Wald tests in [`summary()`](https://rdrr.io/r/base/summary.html)
(individual and joint) to decide which predictors to keep.

## Usage

``` r
# S3 method for class 'ivx'
drop1(
  object,
  scope,
  scale = 0,
  all.cols = TRUE,
  test = c("none", "Chisq", "F"),
  k = 2,
  ...
)

# S3 method for class 'ivx'
add1(
  object,
  scope,
  scale = 0,
  test = c("none", "Chisq", "F"),
  x = NULL,
  k = 2,
  ...
)
```

## Arguments

- object:

  an object of class "ivx".

- scope, scale, all.cols, test, k, x, ...:

  as in [`stats::drop1()`](https://rdrr.io/r/stats/add1.html) and
  [`stats::add1()`](https://rdrr.io/r/stats/add1.html).

## Value

as in [`stats::drop1.lm()`](https://rdrr.io/r/stats/add1.html) /
[`stats::add1.lm()`](https://rdrr.io/r/stats/add1.html), with a
deprecation warning.

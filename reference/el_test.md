# Unified Empirical Likelihood Test for Predictability (Liu, Yang, Cai & Peng)

`el_test` implements the unified empirical likelihood (EL) test of Liu,
Yang, Cai and Peng (2019, Section 2.2) for the predictive regression
with unknown intercept \$\$Y_t = \alpha + \beta_1 \Delta X\_{t-1} +
\beta_2 X\_{t-2} + U_t,\$\$ where the lagged difference of the predictor
is included so that the response can be stationary when the predictor is
not. The intercept is removed by differencing at lag \\m = \lfloor n/2
\rfloor\\ (\\\tilde Y_t = Y\_{t+m} - Y_t\\, \\\tilde X_t = X\_{t+m} -
X_t\\), and the EL function is built from the score equations \\\tilde
Z\_{t1} = \tilde e_t \Delta\tilde X\_{t-1}\\ and \\\tilde Z\_{t2} =
\tilde e_t \tilde X\_{t-2}/\sqrt{1 + \tilde X\_{t-2}^2}\\, \\t = 3,
\dots, m\\, with \\\tilde e_t\\ the model error. The weight on the
second equation makes its sample variance converge whether the predictor
is stationary, nearly integrated or a unit root, so that the profile EL
ratios for \\H_0: \beta_2 = 0\\ (no predictability), \\H_0: \beta_1 =
0\\ and the joint null are \\\chi^2(1)\\, \\\chi^2(1)\\ and
\\\chi^2(2)\\ without knowing the persistence (Theorem 2). Single
predictor only.

## Usage

``` r
el_test(formula, data, na.action)

el_test_fit(y, x)

# S3 method for class 'el_test'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted.

- data:

  n optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in data,
  the variables are taken from environment(formula), typically the
  environment from which lm is called.

- na.action:

  a function which indicates what should happen when the data contain
  NAs. The default is set by the na.action setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

- y:

  response vector.

- x:

  an object of class "el_test".

- digits:

  minimal number of significant digits.

- ...:

  unused.

## Value

an object of class "el_test": a list with the EL ratio statistics `stat`
(named `beta2`, `beta1`, `joint`), their `p.value`, the unconstrained EL
estimates `estimate`, the OLS estimates `ols` of \\(\beta_1, \beta_2)\\
on the differenced data, and `m`.

## References

Liu, X., Yang, B., Cai, Z., & Peng, L. (2019). A unified test for
predictability of asset returns regardless of properties of predicting
variables. Journal of Econometrics, 208(1), 141-159.

Owen, A. B. (2001). Empirical Likelihood. Chapman & Hall.

## Examples

``` r
el_test(Ret ~ DP, data = kms)
#> 
#> Call:
#> el_test(formula = Ret ~ DP, data = kms)
#> 
#> Unified empirical likelihood test (Liu, Yang, Cai & Peng, 2019), m = 516
#> 
#>                  Estimate EL ratio    df Pr(> chi)
#> beta1 (dX[t-1]) -0.133221    4.456 1.000    0.0348
#> beta2 (X[t-2])   0.008196    2.404 1.000    0.1210
#> joint                        7.250 2.000    0.0266
#> 
```

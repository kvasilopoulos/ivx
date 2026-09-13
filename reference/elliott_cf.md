# Control-Function Predictability Test of Elliott (2011)

`elliott_cf` runs the augmented predictive regression of Elliott (2011,
eq. 9), \$\$y_t = \alpha + \beta' x\_{t-1} + \gamma' Z_t + \tilde
u_t,\qquad Z_t = (z_t', z\_{t-1}', \dots, z\_{t-q}')',\$\$ where \\z_t\\
are user-supplied stationary covariates that are contemporaneously
correlated with the innovations of the persistent predictors \\x_t\\ and
of the response. If they absorb that correlation (the "orthogonalising"
condition, Section 5 of the paper), the Wald test of \\\beta = 0\\ has a
standard chi-square limit whatever the persistence of \\x_t\\ (Theorem
2); without them it has the non-standard Elliott-Stock (1994)
distribution (Theorem 1). Unlike
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) and
[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md),
which build the control variable from the data, the covariates here come
from the user

- the paper's example is predicting returns with the dividend-price
  ratio using contemporaneous price-related variables.

## Usage

``` r
elliott_cf(formula, covariates, data, lags = 0, robust = TRUE, na.action)

elliott_cf_fit(y, x, z, lags = 0, robust = TRUE)

# S3 method for class 'elliott_cf'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted.

- covariates:

  a one-sided formula giving the orthogonalising covariates \\z_t\\
  (contemporaneous with the response).

- data:

  n optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in data,
  the variables are taken from environment(formula), typically the
  environment from which lm is called.

- lags:

  number of lags \\q\\ of the covariates to include.

- robust:

  logical; if `TRUE` (default) Eicker-White standard errors are used.

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

  an object of class "elliott_cf".

- z:

  matrix of covariates.

- digits:

  minimal number of significant digits.

- ...:

  unused.

## Value

an object of class "elliott_cf": a list with `coefficients` (on the
predictors), `se`, `tstat`, `Wald`, `df`, `p.value`, `gamma`
(coefficients on the covariates and their lags), `rho_resid` (remaining
innovation correlation per predictor) and the underlying `lm` fit.

## Details

The remaining innovation correlation after augmentation is returned as a
diagnostic: it is the correlation between the regression residuals and
the residuals of an AR(1) of each predictor on the same covariates, and
should be close to zero for the test to be reliable.

## References

Elliott, G. (2011). A control function approach for testing the
usefulness of trending variables in predictive regressions and
econometric models. Journal of Econometrics, 164(1), 79-91.

Elliott, G., & Stock, J. H. (1994). Inference in time series regression
when the order of integration of a regressor is unknown. Econometric
Theory, 10(3-4), 672-700.

## Examples

``` r
# the T-bill rate as covariate for the dividend-price ratio (illustration only)
elliott_cf(Ret ~ DP, ~ TBL, data = kms)
#> 
#> Call:
#> elliott_cf(formula = Ret ~ DP, covariates = ~TBL, data = kms)
#> 
#> Control-function predictive regression (Elliott, 2011), 0 covariate lag(s)
#> 
#>    Estimate Std. Error t value Pr(>|t|)
#> DP 0.005664   0.005099   1.111    0.267
#> (Eicker-White standard errors)
#> 
#> Wald statistic: 1.234 on 1 DF, p-value 0.2667
#> Remaining innovation correlation: DP -0.977 
#> 
```

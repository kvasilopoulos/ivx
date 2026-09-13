# Fitting IVX Models

ivx fits predictive regression models. The method allows standard
chi-square testing for regressors with different degrees of persistence,
from stationary to mildly explosive, and can be used for both short- and
long-horizon predictive regressions.

## Usage

``` r
ivx(
  formula,
  data,
  horizon,
  na.action,
  weights,
  contrasts = NULL,
  offset,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  beta = 0.95,
  cz = 1,
  bandwidth = NULL,
  robust = FALSE,
  lag_y = FALSE,
  ...
)

# S3 method for class 'ivx'
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

- horizon:

  is the horizon (default horizon = 1 corresponds to a short-horizon
  regression). For `horizon > 1` the estimator and modified Wald
  statistic are those of Kostakis et al. (2023), eqs (15) and (23):
  K-period sums of the response and of the lagged predictors, a
  single-lag instrument in the signal matrix and the K-period sum of the
  instrument in the covariance.

- na.action:

  a function which indicates what should happen when the data contain
  NAs. The default is set by the na.action setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector. If non-NULL, weighted least
  squares is used with weights `weights` (that is, minimizing
  `sum(w*e^2)`); otherwise ordinary least squares is used.

- contrasts:

  an optional list. See the `contrasts.arg` of
  [`model.matrix.default`](https://rdrr.io/r/stats/model.matrix.html).

- offset:

  this can be used to specify an a priori known component to be included
  in the linear predictor during fitting. This should be NULL or a
  numeric vector or matrix of extents matching those of the response.
  One or more offset terms can be included in the formula instead or as
  well, and if more than one are specified their sum is used. See
  [model.offset](https://rdrr.io/r/stats/model.extract.html)

- model:

  logical. If `TRUE` the model.frame of the fit is returned.

- x:

  an object of class "ivx", usually, a result of a call to ivx.

- y:

  logical. If `TRUE` the response of the fit is returned.

- beta, cz:

  tuning parameters of the IVX instrument \\z_t = \sum\_{j=0}^{t-1} (1 -
  c_z/n^\beta)^j \Delta x\_{t-j}\\. Defaults (`beta = 0.95`, `cz = 1`)
  follow Kostakis et al. (2015).

- bandwidth:

  Newey-West bandwidth for the long-run covariance estimate. The default
  `NULL` uses \\n^{1/3}\\ as in Kostakis et al. (2015).

- robust:

  logical. If `TRUE` the Eicker-White (heteroskedasticity-robust) form
  of the IVX covariance matrix is used (Demetrescu et al., 2023). Only
  available for `horizon = 1`.

- lag_y:

  logical. If `TRUE` the regression is augmented with the lagged
  dependent variable (column `y_lag`), instrumented by itself, as in
  Demetrescu (2014): this can raise the local power of the IVX test when
  the predictors are highly persistent and endogenous, at no cost
  otherwise. The joint Wald statistic still tests only the predictors.
  Only for `horizon = 1`.

- ...:

  additional arguments to be passed to the low level regression fitting
  functions (see [lm](https://rdrr.io/r/stats/lm.html)).

- digits:

  the number of significant digits to use when printing.

## Value

an object of class "ivx".

## References

Magdalinos, T., & Phillips, P. (2009). Limit Theory for Cointegrated
Systems with Moderately Integrated and Moderately Explosive Regressors.
Econometric Theory, 25(2), 482-526.

Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). Robust
econometric inference for stock return predictability. The Review of
Financial Studies, 28(5), 1506-1553.

Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
(2023). Extensions to IVX methods of inference for return
predictability. Journal of Econometrics, 237(2), 105271.

Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023). Taking
stock of long-horizon predictability tests: Are factor returns
predictable? Journal of Econometrics, 237(2), 105380.

Demetrescu, M. (2014). Enhancing the local power of IVX-based tests in
predictive regressions. Economics Letters, 124(2), 269-273.

## Examples

``` r

# Univariate
ivx(Ret ~ LTY, data = kms)
#> 
#> Call:
#> ivx(formula = Ret ~ LTY, data = kms, horizon = 1)
#> 
#> Coefficients:
#>      LTY  
#> -0.06649  
#> 

# Multivariate
ivx(Ret ~ LTY + TBL, data = kms)
#> 
#> Call:
#> ivx(formula = Ret ~ LTY + TBL, data = kms, horizon = 1)
#> 
#> Coefficients:
#>      LTY       TBL  
#>  0.07624  -0.13497  
#> 

# Longer horizon
ivx(Ret ~ LTY + TBL, data = kms, horizon = 4)
#> 
#> Call:
#> ivx(formula = Ret ~ LTY + TBL, data = kms, horizon = 4)
#> 
#> Coefficients:
#>      LTY       TBL  
#>  0.09322  -0.14164  
#> 

wt <- runif(nrow(kms))
ivx(Ret ~ LTY, data = kms, weights = wt)
#> 
#> Call:
#> ivx(formula = Ret ~ LTY, data = kms, weights = wt, horizon = 1)
#> 
#> Coefficients:
#>     LTY  
#> -0.0705  
#> 

# lag-augmented IVX (Demetrescu, 2014)
ivx(Ret ~ DP, data = kms, lag_y = TRUE)
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, lag_y = TRUE, horizon = 1)
#> 
#> Coefficients:
#>       DP     y_lag  
#> 0.007587  0.094223  
#> 
```

# Fitting IVX-AR Models

ivx_ar implements the Yang et al (2020) new instrumental variable based
Wald statistic (IVX-AR) which accounts for serial correlation and
heteroscedasticity in the error terms of the linear predictive
regression model.

## Usage

``` r
ivx_ar(
  formula,
  data,
  horizon,
  ar = "auto",
  ar_ic = c("bic", "aic", "aicc"),
  ar_max = 5,
  ar_grid = function(x) seq(x - 0.3, x + 0.3, by = 0.02),
  na.action,
  contrasts = NULL,
  offset,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  beta = 0.95,
  cz = 1,
  bandwidth = NULL,
  robust = FALSE,
  ...
)

# S3 method for class 'ivx_ar'
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

- ar:

  Method to include the autoregressive terms. "auto" find the optimal ar
  order by using the information criteria. `ar = 0` reduces to simple
  [`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md).
  `ar > 1` uses a fixed order to estimate the model.

- ar_ic:

  Information criterion to be used in model selection.

- ar_max:

  Maximum ar order of model to fit.

- ar_grid:

  The ar grid sequence of which to iterate.

- na.action:

  a function which indicates what should happen when the data contain
  NAs. The default is set by the na.action setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

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

  an object of class "ivx_ar", usually, a result of a call to ivx_ar.

- y:

  logical. If `TRUE` the response of the fit is returned.

- beta, cz:

  tuning parameters of the IVX instrument \\z_t = \sum\_{j=0}^{t-1} (1 -
  c_z/n^eta)^j \Delta x\_{t-j}\\. Defaults (`beta = 0.95`, `cz = 1`)
  follow Kostakis et al. (2015).

- bandwidth:

  Newey-West bandwidth for the long-run covariance estimate. The default
  `NULL` uses \\n^{1/3}\\ as in Kostakis et al. (2015).

- robust:

  logical. If `TRUE` the Eicker-White (heteroskedasticity-robust) form
  of the IVX covariance matrix is used (Demetrescu et al., 2023). Only
  available for `horizon = 1`.

- ...:

  additional arguments to be passed to the low level regression fitting
  functions (see [lm](https://rdrr.io/r/stats/lm.html)).

- digits:

  the number of significant digits to use when printing.

## References

Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the
Predictability of US Housing Price Index Returns Based on an IVX-AR
Model. Journal of the American Statistical Association, 1-22. DOI:
10.1080/01621459.2019.1686392

## Examples

``` r

ivx_ar(hpi ~ log(res) + cpi, ylpc)
#> 
#> Call:
#> ivx_ar(formula = hpi ~ log(res) + cpi, data = ylpc, horizon = 1)
#> 
#> Lag Selection:
#> Auto (bic) with AR terms q = 4
#> 
#> Coefficients:
#>   log(res)         cpi  
#>  0.0018317  -0.0002078  
#> 

ivx_ar(hpi ~ log(res) + cpi, ylpc, ar_ic = "aic")
#> 
#> Call:
#> ivx_ar(formula = hpi ~ log(res) + cpi, data = ylpc, ar_ic = "aic", 
#>     horizon = 1)
#> 
#> Lag Selection:
#> Auto (aic) with AR terms q = 4
#> 
#> Coefficients:
#>   log(res)         cpi  
#>  0.0018317  -0.0002078  
#> 

ivx_ar(hpi ~ log(res) + cpi, ylpc, ar = 1)
#> 
#> Call:
#> ivx_ar(formula = hpi ~ log(res) + cpi, data = ylpc, ar = 1, horizon = 1)
#> 
#> Lag Selection:
#> Fixed with AR terms q = 1
#> 
#> Coefficients:
#>   log(res)         cpi  
#>  0.0010780  -0.0001155  
#> 
```

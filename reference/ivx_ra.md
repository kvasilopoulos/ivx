# Fitting Residual-Augmented IVX Models

`ivx_ra` implements the residual-augmented (bias-reduced) IVX estimator
of Demetrescu and Rodrigues (2022). An autoregression of order `p` is
fitted to the predictors, the predictive regression is augmented with
its residuals (in the spirit of Amihud and Hurvich, 2004), and the slope
on the lagged predictors is estimated by IVX. Inference uses the
heteroskedasticity-robust standard errors of the paper (eq. 9 and 14),
which are valid whether the predictors are stationary or
near-integrated.

## Usage

``` r
ivx_ra(
  formula,
  data,
  ar = "auto",
  ar_ic = c("aic", "bic"),
  ar_max = 5,
  horizon = 1,
  beta = 0.95,
  cz = 1,
  na.action,
  contrasts = NULL,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  ...
)

# S3 method for class 'ivx_ra'
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

- ar:

  order of the autoregression fitted to the predictors: `"auto"` selects
  it by `ar_ic` in levels (as recommended in the paper), or a positive
  integer for a fixed order.

- ar_ic:

  information criterion for `ar = "auto"`.

- ar_max:

  maximum order considered when `ar = "auto"`.

- horizon:

  forecast horizon \\h\\; see Details.

- beta, cz:

  tuning parameters of the IVX instrument \\z_t = \sum\_{j=0}^{t-1} (1 -
  c_z/n^eta)^j \Delta x\_{t-j}\\. Defaults (`beta = 0.95`, `cz = 1`)
  follow Kostakis et al. (2015).

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

- model:

  logical. If `TRUE` the model.frame of the fit is returned.

- x:

  an object of class "ivx_ra".

- y:

  logical. If `TRUE` the response of the fit is returned.

- ...:

  additional arguments to be passed to the low level regression fitting
  functions (see [lm](https://rdrr.io/r/stats/lm.html)).

- digits:

  the number of significant digits to use when printing.

## Value

an object of class `c("ivx_ra", "ivx")`; the usual `ivx` methods
(`summary`, `vcov`, `delta`, ...) apply. Additional components: `gamma`
(coefficients on the augmentation residuals) and `ar_order`.

## Details

For `horizon > 1` the estimator is the transformed-regression
residual-augmented IVX of Demetrescu, Rodrigues and Taylor (2023), eqs
(4.9), (4.11) and (5.5)-(5.7): the single-period response is regressed
on the \\h\\-period transformed instrument \\z_t^{trf,(h)} =
\sum\_{i=\max(1,t-h+1)}^{\min(t,T-h)} z_i\\ (eq. 4.4), which accounts
for the overlap of the long-horizon regression without HAC estimation.
At `horizon = 1` it coincides with the short-horizon estimator. The
coefficients estimate the \\h\\-period slope \\eta_h\\; fitted values
and residuals are those of the transformed (non-overlapping) regression,
and the Kostakis et al. (2015) intercept correction is applied only at
`horizon = 1`.

The autoregression of the predictors is fitted without an intercept and
its residuals are demeaned before augmentation, which is the paper's
preferred \\\tilde t\_{ivx}^{\mu_0}\\ statistic (Sections 4-5); the
standard errors include the finite-sample correction of Kostakis et al.
(2015) as in the paper's simulations.

## References

Demetrescu, M., & Rodrigues, P. M. M. (2022). Residual-augmented IVX
predictive regression. Journal of Econometrics, 227(2), 429-460.

Demetrescu, M., Rodrigues, P. M. M., & Taylor, A. M. R. (2023).
Transformed regression-based long-horizon predictability tests. Journal
of Econometrics, 237(2), 105316.

Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
reduced-bias estimation method. Journal of Financial and Quantitative
Analysis, 39(4), 813-841.

## Examples

``` r
ivx_ra(Ret ~ DP, data = kms)
#> 
#> Call:
#> ivx_ra(formula = Ret ~ DP, data = kms)
#> 
#> Residual-augmented IVX, AR order p = 5 (aic)
#> 
#> Coefficients:
#>        DP  
#> -0.002383  
#> 

summary(ivx_ra(Ret ~ DP + TBL, data = kms, ar = 2))
#> 
#> Call:
#> ivx_ra(formula = Ret ~ DP + TBL, data = kms, ar = 2)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP  -0.002924   0.005016  -0.583    0.340      0.56
#> TBL -0.051281   0.054856  -0.935    0.874      0.35
#> (Eicker-White standard errors)
#> 
#> Joint Wald statistic:  1.654 on 2 DF, p-value 0.4374
#> Multiple R-squared:  0.02494,    Adjusted R-squared:  0.02304
#> 

# long horizon (Demetrescu, Rodrigues & Taylor, 2023)
ivx_ra(Ret ~ DP, data = kms, horizon = 12)
#> 
#> Call:
#> ivx_ra(formula = Ret ~ DP, data = kms, horizon = 12)
#> 
#> Residual-augmented IVX, AR order p = 5 (aic)
#> 
#> Coefficients:
#>       DP  
#> -0.01775  
#> 
```

# Fitting Systems of IVX Predictive Regressions

`ivx_sys` estimates a system of predictive regressions \\y_t = \mu + A
x\_{t-1} + \varepsilon_t\\ with an \\m\\-vector response and \\r\\
predictors of arbitrary persistence by IVX, and tests linear
restrictions on \\A\\ with the IVX-Wald statistic of Kostakis,
Magdalinos and Stamatogiannis (2023), eqs (15) and (23), whose
covariance has the Kronecker form \\(Z'X)^{-1} \otimes I_m\\ around
\\Z(K)'Z(K) \otimes \hat\Sigma - n \bar z \bar z' \otimes
\hat\Sigma\_{FM}\\. Long horizons are handled as in
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md). For \\m
= 1\\ the results coincide with
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md).

## Usage

``` r
ivx_sys(
  formula,
  data,
  horizon,
  na.action,
  contrasts = NULL,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  beta = 0.95,
  cz = 1,
  bandwidth = NULL,
  ...
)

# S3 method for class 'ivx_sys'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'ivx_sys'
summary(object, ...)
```

## Arguments

- formula:

  a formula whose left-hand side is a matrix, e.g.
  `cbind(y1, y2) ~ x1 + x2`.

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

- contrasts:

  an optional list. See the `contrasts.arg` of
  [`model.matrix.default`](https://rdrr.io/r/stats/model.matrix.html).

- model:

  logical. If `TRUE` the model.frame of the fit is returned.

- x:

  an object of class "ivx_sys".

- y:

  logical. If `TRUE` the response of the fit is returned.

- beta, cz:

  tuning parameters of the IVX instrument \\z_t = \sum\_{j=0}^{t-1} (1 -
  c_z/n^eta)^j \Delta x\_{t-j}\\. Defaults (`beta = 0.95`, `cz = 1`)
  follow Kostakis et al. (2015).

- bandwidth:

  Newey-West bandwidth for the long-run covariance estimate. The default
  `NULL` uses \\n^{1/3}\\ as in Kostakis et al. (2015).

- ...:

  additional arguments to be passed to the low level regression fitting
  functions (see [lm](https://rdrr.io/r/stats/lm.html)).

- digits:

  the number of significant digits to use when printing.

- object:

  an object of class "ivx_sys".

## Value

an object of class "ivx_sys" with the coefficient matrix `coefficients`
(responses in rows, predictors in columns), matching `se`, `tstat` and
`Wald_Ind` matrices, the joint Wald statistic `Wald_Joint`
(\\\chi^2(mr)\\), the per-equation Wald statistics `Wald_Eq`
(\\\chi^2(r)\\), and `vcov` of `vec(coefficients)` (column-major, names
`response:predictor`).

## References

Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023). Taking
stock of long-horizon predictability tests: Are factor returns
predictable? Journal of Econometrics, 237(2), 105380.

## Examples

``` r
ivx_sys(cbind(Ret, DE) ~ DP + TBL, data = kms)
#> 
#> Call:
#> ivx_sys(formula = cbind(Ret, DE) ~ DP + TBL, data = kms, horizon = 1)
#> 
#> Coefficients (responses in rows):
#>      DP         TBL      
#> Ret   0.006145  -0.080717
#> DE    0.311083  -4.014921
#> 

summary(ivx_sys(cbind(Ret, DE) ~ DP + TBL, data = kms, horizon = 4))
#> 
#> Call:
#> ivx_sys(formula = cbind(Ret, DE) ~ DP + TBL, data = kms, horizon = 4)
#> 
#> Coefficients:
#>          Estimate Std. Error t value Wald Ind Pr(> chi)    
#> Ret:DP   0.006579   0.004601   1.430    2.045     0.153    
#> DE:DP    0.316547   0.019563  16.181  261.819    <2e-16 ***
#> Ret:TBL -0.073549   0.058238  -1.263    1.595     0.207    
#> DE:TBL  -3.979948   0.279933 -14.217  202.137    <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Equation Wald statistics on 2 DF:
#>      Wald p-value
#> Ret 3.527  0.1715
#> DE  488.3  <2e-16
#> 
#> Joint Wald statistic:  493.7 on 4 DF, p-value < 2.2e-16 
#> 
```

# Fitting IVX Quantile Predictive Regressions

`ivx_qr` implements the IVX-QR predictability test of Lee (2016): the
\\\tau\\-quantile of the response is regressed on the IVX-filtered
predictors (Section 3.3, Proposition 3.2), giving a test of \\H_0:
\beta\_\tau = 0\\ with a standard chi-square limit whatever the
persistence of the predictors. Requires the quantreg package.

## Usage

``` r
ivx_qr(
  formula,
  data,
  tau = 0.5,
  beta = 0.95,
  cz = 5,
  na.action,
  contrasts = NULL,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  ...
)

# S3 method for class 'ivx_qr'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'ivx_qr'
summary(object, ...)
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

- tau:

  quantile level(s) in (0, 1). A vector fits one model per level.

- beta, cz:

  tuning parameters of the IVX instrument \\z_t = \sum\_{j=0}^{t-1} (1 -
  c_z/n^\beta)^j \Delta x\_{t-j}\\. Defaults (`beta = 0.95`, `cz = 1`)
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

  an object of class "ivx_qr".

- y:

  logical. If `TRUE` the response of the fit is returned.

- ...:

  further arguments passed to
  [`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html).

- digits:

  the number of significant digits to use when printing.

- object:

  an object of class "ivx_qr".

## Value

For a single `tau`, an object of class `c("ivx_qr", "ivx")` (so
[`summary()`](https://rdrr.io/r/base/summary.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) etc. apply) with
components `tau`, `rho_tau`, `sparsity` (\\\hat f_u(0)\\) and `rq` (the
underlying [`quantreg::rq`](https://rdrr.io/pkg/quantreg/man/rq.html)
fit). For several `tau`, a list of such objects named by `tau`.

## Details

The instrument is \\z_t = (1 - c_z/n^\beta) z\_{t-1} + \Delta x_t\\. Lee
(2016) normalises \\c_z = 5\\ and picks \\\beta\\ from a look-up table
indexed by the estimated QR endogeneity \\\hat\rho(\tau) =
-\mathrm{corr}(1\\\hat u_t \< 0\\, \hat u\_{x,t})\\, which is returned
as `rho_tau` so the rule can be applied by the user; the default
`beta = 0.95` follows Kostakis et al. (2015). The sparsity \\f_u(0)\\ is
estimated by a Gaussian kernel with Silverman's bandwidth (footnote 4 of
the paper).

## References

Lee, J. H. (2016). Predictive quantile regression with persistent
covariates: IVX-QR approach. Journal of Econometrics, 192(1), 105-118.

## Examples

``` r
if (requireNamespace("quantreg", quietly = TRUE)) {
  summary(ivx_qr(Ret ~ DP, data = kms, tau = 0.5))
  ivx_qr(Ret ~ DP + TBL, data = kms, tau = c(0.1, 0.5, 0.9))
}
#> $`0.1`
#> 
#> Call:
#> ivx_qr(formula = Ret ~ DP + TBL, data = kms, tau = 0.1)
#> 
#> IVX-QR at tau = 0.1
#> 
#> Coefficients:
#>       DP       TBL  
#> -0.02831   0.18938  
#> 
#> 
#> $`0.5`
#> 
#> Call:
#> ivx_qr(formula = Ret ~ DP + TBL, data = kms, tau = 0.5)
#> 
#> IVX-QR at tau = 0.5
#> 
#> Coefficients:
#>        DP        TBL  
#>  0.008204  -0.171321  
#> 
#> 
#> $`0.9`
#> 
#> Call:
#> ivx_qr(formula = Ret ~ DP + TBL, data = kms, tau = 0.9)
#> 
#> IVX-QR at tau = 0.9
#> 
#> Coefficients:
#>       DP       TBL  
#>  0.02674  -0.43021  
#> 
#> 
```

# Wild Bootstrap Inference for IVX Models

Computes bootstrap p-values for the IVX Wald and t statistics using the
residual wild bootstrap (RWB) or the fixed regressor wild bootstrap
(FRWB) of Demetrescu et al. (2023). The null hypothesis of no
predictability is imposed on the bootstrap samples. RWB rebuilds the
regressor from an AR fit and its residuals (multiplied by the same wild
multiplier as the predictive-regression residuals, so the innovation
correlation is preserved); FRWB keeps the regressors and instruments
fixed and only resamples the response.

## Usage

``` r
ivx_boot(
  object,
  B = 999,
  type = c("rwb", "frwb"),
  ar_max = 5,
  dist = c("rademacher", "normal"),
  seed = NULL,
  cores = 1L
)

# S3 method for class 'ivx_boot'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  an object of class "ivx" (not "ivx_ar"), fitted without weights.

- B:

  number of bootstrap replications.

- type:

  bootstrap scheme: `"rwb"` (residual wild bootstrap, recommended for
  strongly persistent regressors) or `"frwb"` (fixed regressor wild
  bootstrap).

- ar_max:

  maximum lag order of the (vector) autoregression fitted to the
  regressors by the RWB scheme; the order is selected by BIC (Remark
  24).

- dist:

  distribution of the wild multipliers.

- seed:

  optional integer seed. With `cores > 1` the L'Ecuyer-CMRG streams of
  the parallel package are used, so results are reproducible for a given
  `seed` and `cores` but differ from the serial run.

- cores:

  number of CPU cores. Uses forking on Unix and a PSOCK cluster on
  Windows (the package must be installed for the workers to load it).

- x:

  an object of class "ivx_boot".

- digits:

  minimal number of significant digits.

- ...:

  further arguments passed to
  [`printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html).

## Value

an object of class "ivx_boot": a list with the observed statistics, the
bootstrap distributions (`boot`), and bootstrap p-values (`p.value`).
`p.value$tstat` has one column per alternative: two-sided, `beta < 0`
and `beta > 0`.

## References

Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
(2023). Extensions to IVX methods of inference for return
predictability. Journal of Econometrics, 237(2), 105271.

## Examples

``` r
mod <- ivx(Ret ~ DP + TBL, data = kms)
ivx_boot(mod, B = 199, seed = 1)
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 1)
#> 
#> Residual wild bootstrap, B = 199
#> 
#> Coefficients (bootstrap p-values):
#>      Estimate t value Wald Ind Pr(> chi) Pr(t < 0) Pr(t > 0)
#> DP   0.006145   1.349    1.819    0.3668    0.6935     0.307
#> TBL -0.080717  -1.399    1.957    0.2060    0.0804     0.920
#> 
#> Joint Wald statistic: 3.644, bootstrap p-value 0.3266
#> 
```

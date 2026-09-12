# Subsample IVX Tests for Episodic Predictability

Tests for "pockets" of predictability using the suprema of sequences of
subsample IVX statistics (Demetrescu et al., 2023, Section 3.2;
Demetrescu et al., 2022). For each window the IVX statistic is computed
from the window's observations with the full-sample instrument (eqs
15-17); the test statistics are the maximum (right-tailed), minimum
(left-tailed) and maximum squared (two-sided) t-ratio over the sequence
for a single predictor, and the maximum Wald statistic for several
predictors (Remark 11). P-values are obtained by wild bootstrap
(Algorithms 1-2), which the paper shows to be asymptotically valid for
these sup-functionals.

## Usage

``` r
ivx_episodic(
  object,
  scheme = c("rolling", "forward", "backward"),
  window = 0.2,
  robust = FALSE,
  B = 999,
  type = c("frwb", "rwb"),
  ar_max = 5,
  dist = c("rademacher", "normal"),
  seed = NULL,
  cores = 1L
)

# S3 method for class 'ivx_episodic'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  an object of class "ivx" fitted with `horizon = 1`.

- scheme:

  `"rolling"` windows of fixed width, `"forward"` recursive windows
  starting at the first observation, or `"backward"` recursive windows
  ending at the last observation.

- window:

  fraction of the sample: the window width for `"rolling"`, the warm-in
  fraction \\\tau_L\\ for `"forward"`, and the latest start \\\tau_U\\
  for `"backward"`.

- robust:

  logical; use Eicker-White standard errors in the subsample statistics.

- B:

  number of bootstrap replications.

- type:

  bootstrap scheme, see
  [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md).
  The default fixed regressor wild bootstrap is the scheme used by
  Demetrescu et al. (2022).

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

  an object of class "ivx_episodic".

- digits:

  minimal number of significant digits.

- ...:

  unused.

## Value

an object of class "ivx_episodic": the observed sequence (`sequence`,
one row per window with its start/end index and statistics), the sup
statistics (`statistic`) and their bootstrap p-values (`p.value`), plus
the bootstrap draws (`boot`).

## References

Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
(2022). Testing for episodic predictability in stock returns. Journal of
Econometrics, 227(1), 85-113.

Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
(2023). Extensions to IVX methods of inference for return
predictability. Journal of Econometrics, 237(2), 105271.

## Examples

``` r
mod <- ivx(Ret ~ DP, data = kms)
ivx_episodic(mod, scheme = "rolling", window = 0.2, B = 99, seed = 1)
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, horizon = 1)
#> 
#> Subsample IVX tests, rolling scheme (window = 0.2), 827 windows
#> Fixed regressor wild bootstrap, B = 99
#> 
#>                         statistic bootstrap p
#> sup t   (H1: beta > 0)      2.921       0.202
#> inf t   (H1: beta < 0)    -0.4166       1.000
#> sup t^2 (H1: beta != 0)     8.533       0.303
#> 
```

# IV Predictability Tests with Combined Instruments

`ivx_iv` implements the instrumental-variable tests of Breitung and
Demetrescu (2015): the predictive regression is estimated by 2SLS with
Eicker-White standard errors (their eq. 12) using instruments that are
less persistent than the predictor. Two families are available. Type-I
instruments are transformations of the predictor itself: the fractional
difference \\\Delta\_+^{d} x\_{t-1}\\ (`"frac"`) and the long difference
\\x\_{t-1} - x\_{t-1-k_T}\\ (`"diff"`), which keep power when the
predictor is stationary. Type-II instruments are deterministic and
correlate with a near-integrated predictor only: the sine function
\\\sin(\pi t/T)\\ (`"sin"`). Their 2SLS combination (`"comb"`, the
paper's `IVcomb` and the authors' recommendation) is asymptotically
dominated by whichever instrument is informative, so the squared t-ratio
(and the Wald statistic with several predictors) is chi-square whatever
the persistence of the predictor (Theorems 3 and 7).

## Usage

``` r
ivx_iv(
  formula,
  data,
  instruments = c("comb", "sin", "frac", "diff"),
  d = 0.5,
  kappa = 0.2,
  eta = 0.85,
  na.action,
  contrasts = NULL,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  ...
)

# S3 method for class 'ivx_iv'
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

- instruments:

  instrument set; see Details.

- d:

  order of the fractional difference for `"frac"`, in (0, 1/2\]; the
  paper uses 1/2.

- kappa, eta:

  the long-difference lag is \\k_T = \lfloor \kappa T^\eta \rfloor\\
  (paper: 0.2 and 0.85), truncated to \\t - 1\\.

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

  an object of class "ivx_iv".

- y:

  logical. If `TRUE` the response of the fit is returned.

- ...:

  additional arguments to be passed to the low level regression fitting
  functions (see [lm](https://rdrr.io/r/stats/lm.html)).

- digits:

  the number of significant digits to use when printing.

## Value

an object of class `c("ivx_iv", "ivx")` so the `ivx` methods apply;
`instruments` holds the instrument matrix aligned with the regressors.

## Details

With \\K\\ predictors each type-I instrument is built per predictor and
the sine instruments use frequencies \\\sin(k\pi t/T)\\, \\k = 1, \dots,
K\\, so that the instrument vector stays linearly independent
(Assumption 5). The IVX instrument of Kostakis et al. (2015) is the
paper's "mild integration" type-I case and is available through
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md).

## References

Breitung, J., & Demetrescu, M. (2015). Instrumental variable and
variable addition based inference in predictive regressions. Journal of
Econometrics, 187(1), 358-375.

## Examples

``` r
ivx_iv(Ret ~ DP, data = kms)
#> 
#> Call:
#> ivx_iv(formula = Ret ~ DP, data = kms)
#> 
#> 2SLS with instruments: sin1, frac_DP
#> 
#> Coefficients:
#>     DP  
#> 0.0101  
#> 
summary(ivx_iv(Ret ~ DP + TBL, data = kms, instruments = "frac"))
#> 
#> Call:
#> ivx_iv(formula = Ret ~ DP + TBL, data = kms, instruments = "frac")
#> 
#> Coefficients:
#>       Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP  -0.0455575  0.2025258  -0.225    0.051     0.822
#> TBL -0.0009177  0.8732213  -0.001    0.000     0.999
#> (Eicker-White standard errors)
#> 
#> Joint Wald statistic:  3.117 on 2 DF, p-value 0.2104
#> Multiple R-squared:  0.8671, Adjusted R-squared:  0.8667
#> 
```

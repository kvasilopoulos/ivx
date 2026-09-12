# Hybrid t-Test for Return Predictability (Harvey, Leybourne & Taylor)

`hlt_test` implements the double-switching hybrid procedure \\T\_{hyb}\\
of Harvey, Leybourne and Taylor (2021, Section 3.3) for a single
predictor. Two regression t-ratios are used: the standard OLS t-ratio
\\T\\ (eq. 5) and the variant \\\tilde T\\ in which the predictor is
quasi-GLS demeaned with \\\bar c = 7\\ (eq. 7). Under weak persistence
the standard t-ratio is compared with normal critical values; under
strong persistence the limiting null distributions depend on the
local-to-unity parameter and on the innovation correlation
\\\rho\_{xy}\\, and the tests are run with the paper's asymptotically
conservative critical values (maximised over \\c\\) obtained from the
response surfaces in its Table 1.

## Usage

``` r
hlt_test(
  formula,
  data,
  alternative = c("greater", "less"),
  level = 0.05,
  lag_max = NULL,
  na.action
)

hlt_test_fit(y, x, alternative = "greater", level = 0.05, lag_max = NULL)

# S3 method for class 'hlt_test'
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

- alternative:

  direction of the one-sided test.

- level:

  significance level; one of 0.1, 0.05, 0.025, 0.01 (the levels for
  which response surfaces are available).

- lag_max:

  maximum ADF lag order; the default is the paper's \\\lfloor 12
  (T/100)^{1/4} \rfloor\\.

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

  an object of class "hlt_test".

- digits:

  minimal number of significant digits.

- ...:

  unused.

## Value

an object of class "hlt_test": a list with the selected `test` (`"T_N"`,
`"T_con"` or `"T~_con"`), its `statistic`, `cv` and `reject` indicator,
plus `t`, `t_gls` (both t-ratios), `adf`, `adf_lag`, `rho_xy` and
`estimate` (the OLS slope).

## Details

The procedure is: (1) if the ADF normalised-bias statistic
\\T\hat\rho/(1 - \sum_i \hat\psi_i)\\ from an ADF regression with lag
length chosen by the MBIC of Ng and Perron (2001) is below
\\-4\sqrt{T}\\, the predictor is treated as weakly persistent and the
standard test \\T_N\\ is used; (2) otherwise, for an upper-tail test,
\\T\\ with critical value \\cv(\hat\rho\_{xy})\\ if \\\hat\rho\_{xy} \>
-0.1\\ and \\\tilde T\\ with \\\tilde{cv}(\hat\rho\_{xy})\\ if
\\\hat\rho\_{xy} \< -0.1\\ (mirrored for lower-tail tests), where
\\\hat\rho\_{xy}\\ is the correlation of the ADF residuals with the
predictive-regression residuals.

## References

Harvey, D. I., Leybourne, S. J., & Taylor, A. M. R. (2021). Simple tests
for stock return predictability with good size and power properties.
Journal of Econometrics, 224(1), 198-214.

Ng, S., & Perron, P. (2001). Lag length selection and the construction
of unit root tests with good size and power. Econometrica, 69(6),
1519-1554.

## Examples

``` r
hlt_test(Ret ~ DP, data = kms)
#> 
#> Call:
#> hlt_test(formula = Ret ~ DP, data = kms)
#> 
#> Hybrid predictability test (Harvey, Leybourne & Taylor, 2021)
#> 
#> Selected test: T~_con (quasi-GLS-demeaned t, conservative critical value)
#> statistic = 1.298, 5% critical value = 1.944 (alternative: beta > 0): do not reject
#> slope = 0.006172, t =  1.63, t (quasi-GLS) = 1.298, ADF = -5.97 (p = 21, cutoff -128.6), rho_xy = -0.9485
#> 
hlt_test(Ret ~ TBL, data = kms, alternative = "less", level = 0.1)
#> 
#> Call:
#> hlt_test(formula = Ret ~ TBL, data = kms, alternative = "less", 
#>     level = 0.1)
#> 
#> Hybrid predictability test (Harvey, Leybourne & Taylor, 2021)
#> 
#> Selected test: T_con (OLS-demeaned t, conservative critical value)
#> statistic = -1.403, 10% critical value = -1.306 (alternative: beta < 0): reject
#> slope = -0.07836, t = -1.403, t (quasi-GLS) = -1.384, ADF = -6.707 (p = 6, cutoff -128.6), rho_xy = -0.05541
#> 
```

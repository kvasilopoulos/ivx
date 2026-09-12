# Augmented Regression Method (Amihud, Hurvich & Wang)

`arm` implements the multipredictor augmented regression method (mARM)
of Amihud, Hurvich and Wang (2009), a reduced-bias OLS alternative to
IVX for stationary but persistent predictors. A VAR(1) is fitted to the
predictors, its coefficient matrix is bias-corrected with the Nicholls
and Pope (1988) expansion (iterated), and the predictive regression is
augmented with the corrected VAR residuals, which removes the Stambaugh
(1999) bias from the slopes. Standard errors and the joint Wald
statistic use the covariance estimator of the paper (eqs 7-8), which
adds the estimation uncertainty of the VAR coefficients to the
augmented-regression OLS variance.

## Usage

``` r
arm(
  formula,
  data,
  iter = 10,
  na.action,
  contrasts = NULL,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  ...
)

# S3 method for class 'arm'
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

- iter:

  maximum number of bias-correction iterations (`K = 10` in the paper);
  iteration stops earlier if the corrected VAR becomes non-stationary.

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

  an object of class "arm".

- y:

  logical. If `TRUE` the response of the fit is returned.

- ...:

  additional arguments to be passed to the low level regression fitting
  functions (see [lm](https://rdrr.io/r/stats/lm.html)).

- digits:

  the number of significant digits to use when printing.

## Value

an object of class `c("arm", "ivx")`, so the `ivx` methods apply.
Additional components: `phi` (coefficients on the augmentation
residuals), `Phi` (bias-corrected VAR(1) coefficient matrix, equations
by row) and `Phi_ols`.

## Details

Unlike IVX the method assumes stationary predictors (all eigenvalues of
the VAR coefficient matrix inside the unit circle) and Gaussian
innovations; it is the natural benchmark for the "control function"
approach of Elliott (2011). Short horizon only.

## References

Amihud, Y., Hurvich, C. M., & Wang, Y. (2009). Multiple-predictor
regressions: Hypothesis testing. The Review of Financial Studies, 22(1),
413-434.

Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
reduced-bias estimation method. Journal of Financial and Quantitative
Analysis, 39(4), 813-841.

Nicholls, D. F., & Pope, A. L. (1988). Bias in the estimation of
multivariate autoregressions. Australian Journal of Statistics, 30A,
296-309.

## Examples

``` r
arm(Ret ~ DP, data = kms)
#> 
#> Call:
#> arm(formula = Ret ~ DP, data = kms)
#> 
#> Augmented regression method (reduced-bias OLS)
#> 
#> Coefficients:
#>       DP  
#> 0.002463  
#> 

summary(arm(Ret ~ DP + TBL, data = kms))
#> 
#> Call:
#> arm(formula = Ret ~ DP + TBL, data = kms)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP   0.002050   0.003797   0.540    0.292     0.589
#> TBL -0.046509   0.055999  -0.831    0.690     0.406
#> 
#> Joint Wald statistic:  1.068 on 2 DF, p-value 0.5862
#> Multiple R-squared:  0.02415,    Adjusted R-squared:  0.0194
#> 
```

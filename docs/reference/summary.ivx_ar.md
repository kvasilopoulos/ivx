# Summarizing IVX-AR Model Fits

summary method for class "ivx".

## Usage

``` r
# S3 method for class 'ivx_ar'
summary(object, ...)

# S3 method for class 'summary.ivx_ar'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object:

  object of class "ivx_ar", usually, a result of a call to ivx_ar.

- ...:

  further arguments passed to or from other methods.

- x:

  an object of class `"summary.lm"`, usually, a result of a call to
  `summary.lm`.

- digits:

  the number of significant digits to use when printing.

- signif.stars:

  logical. If `TRUE`, ‘significance stars’ are printed for each
  coefficient.

## Examples

``` r
mod <- ivx_ar(Ret ~ LTY, data = kms)

summary(mod)
#> 
#> Call:
#> ivx_ar(formula = Ret ~ LTY, data = kms, horizon = 1)
#> 
#> Auto () with AR terms q = 5
#> 
#> Coefficients:
#>     Estimate Std. Error t value Wald Ind Pr(> chi)
#> LTY -0.06808    0.07217  -0.943    0.890     0.346
#> 
#> Joint Wald statistic:  0.8897 on 1 DF, p-value 0.3456
#> Multiple R-squared:  0.0009652,  Adjusted R-squared:  -8.522e-06
#> Wald AR statistic: 9.517 on 5 DF, p-value 0.09013
#> 
```

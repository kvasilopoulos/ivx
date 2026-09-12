# Summarizing IVX Model Fits

summary method for class "ivx".

## Usage

``` r
# S3 method for class 'ivx'
summary(object, ...)

# S3 method for class 'summary.ivx'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object:

  object of class "ivx", usually, a result of a call to ivx.

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
mod <- ivx(Ret ~ LTY, data = monthly)

summary(mod)
#> 
#> Call:
#> ivx(formula = Ret ~ LTY, data = monthly, horizon = 1)
#> 
#> Coefficients:
#>     Estimate Std. Error t value Wald Ind Pr(> chi)
#> LTY -0.06649    0.06447  -1.031    1.064     0.302
#> 
#> Joint Wald statistic:  1.064 on 1 DF, p-value 0.3024
#> Multiple R-squared:  0.001133,   Adjusted R-squared:  0.0001645
#> 
```

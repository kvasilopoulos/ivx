# Calculate the delta coefficient

Computes the long-run correlation coefficient between the residuals of
the predictive regression and the autoregressive model for the
regressor.

## Usage

``` r
delta(object)
```

## Arguments

- object:

  on object of class "ivx"

## Value

A vector of the estimated correlation coefficients. This should have row
and column names corresponding to the parameter names given by the coef
method.

## Examples

``` r
mod <- ivx(Ret ~ LTY, data = monthly)

delta(mod)
#> [1] -0.1081335
```

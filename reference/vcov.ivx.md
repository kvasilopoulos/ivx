# Calculate Variance-Covariance Matrix for a Fitted Model Object

Calculate Variance-Covariance Matrix for a Fitted Model Object

## Usage

``` r
# S3 method for class 'ivx'
vcov(object, complete = TRUE, ...)

# S3 method for class 'summary.ivx'
vcov(object, complete = TRUE, ...)
```

## Arguments

- object:

  a fitted ivx and summary.ivx object.

- complete:

  logical indicating if the full variance-covariance matrix should be
  returned. When complete = TRUE, vcov() is compatible with coef().

- ...:

  additional arguments for method functions.

## Value

A matrix of the estimated covariances between the parameter estimates of
the model. This should have row and column names corresponding to the
parameter names given by the coef method.

## Examples

``` r
mod <- ivx(Ret ~ LTY, data = monthly)

vcov(mod)
#>             LTY
#> LTY 0.004156707
```

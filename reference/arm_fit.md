# Fitter Function for the Augmented Regression Method

Basic function called by `arm`. Should only be used directly by
experienced users.

## Usage

``` r
arm_fit(y, x, iter = 10, ...)
```

## Arguments

- y:

  vector of observations of length `n`, or a matrix with `n` rows.

- x:

  design matrix of dimension `n * p`.

- iter:

  maximum number of bias-correction iterations (`K = 10` in the paper);
  iteration stops earlier if the corrected VAR becomes non-stationary.

- ...:

  currently unused.

## Examples

``` r
arm_fit(kms$Ret, as.matrix(kms$DP))$tstat
#>        x1 
#> 0.6507685 
```

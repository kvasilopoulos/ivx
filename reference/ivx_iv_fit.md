# Fitter Function for IV Predictability Tests

Basic function called by `ivx_iv`. Should only be used directly by
experienced users.

## Usage

``` r
ivx_iv_fit(y, x, instruments = "comb", d = 0.5, kappa = 0.2, eta = 0.85, ...)
```

## Arguments

- y:

  vector of observations of length `n`, or a matrix with `n` rows.

- x:

  design matrix of dimension `n * p`.

- instruments:

  instrument set; see Details.

- d:

  order of the fractional difference for `"frac"`, in (0, 1/2\]; the
  paper uses 1/2.

- kappa, eta:

  the long-difference lag is \\k_T = \lfloor \kappa T^\eta \rfloor\\
  (paper: 0.2 and 0.85), truncated to \\t - 1\\.

- ...:

  currently unused.

## Examples

``` r
ivx_iv_fit(kms$Ret, as.matrix(kms$DP))$tstat
#>       x1 
#> 0.490082 
```

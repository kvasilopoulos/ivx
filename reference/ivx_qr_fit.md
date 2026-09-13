# Fitter Function for IVX-QR Models

Basic function called by `ivx_qr`. Should only be used directly by
experienced users.

## Usage

``` r
ivx_qr_fit(y, x, tau = 0.5, beta = 0.95, cz = 5, ...)
```

## Arguments

- y:

  vector of observations of length `n`, or a matrix with `n` rows.

- x:

  design matrix of dimension `n * p`.

- tau:

  quantile level(s) in (0, 1). A vector fits one model per level.

- beta, cz:

  tuning parameters of the IVX instrument \\z_t = \sum\_{j=0}^{t-1} (1 -
  c_z/n^\beta)^j \Delta x\_{t-j}\\. Defaults (`beta = 0.95`, `cz = 1`)
  follow Kostakis et al. (2015).

- ...:

  currently disregarded.

## Examples

``` r
if (requireNamespace("quantreg", quietly = TRUE)) {
  ivx_qr_fit(kms$Ret, as.matrix(kms$DP), tau = 0.5)$Wald_Joint
}
#> [1] 0.4190654
```

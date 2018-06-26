
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ivx

The goal of ivx is to …

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kvasilopoulos/ivx")
```

## Example

``` r
library(ivx)
```

This is a basic example, lets load the data first:

``` r
# Monthly data from Kostakis et al (2014)
X <- monthly[, 2, drop = F]
Y <- monthly[, 13, drop = F]
reg <- data.frame(X = X$`D/E`, Y = Y$LOG_EXCESS_VW)
```

And then do the estimation:

``` r
est <- ivx(Y ~ X, data = reg)
est
#> 
#> Call:
#> ivx(formula = Y ~ X, data = reg)
#> 
#> Coefficients:
#> [1]  -0.003287
```

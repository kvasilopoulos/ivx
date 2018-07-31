
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

Load ivx and dplyr for data managment

``` r
library(ivx)
library(dplyr)
```

This is a basic example, lets load the data first:

``` r
# Monthly data from Kostakis et al (2014)
datam %>%
  names()
#>  [1] "Date" "D_E"  "LTY"  "D_Y"  "D_P"  "TBL"  "E_P"  "B_M"  "INF"  "DFY" 
#> [11] "NTIS" "TMS"  "Ret"

univariate <- datam %>%
  select(Ret, D_P)

multivariate <- datam %>%
  select(Ret, D_P, TBL)
```

And then do the univariate estimation:

``` r
ivx(Ret ~ D_P, data = univariate) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ D_P, data = univariate, horizon = 1)
#> 
#> Coefficients:
#>      Estimate  Wald Ind Pr(> chi)     delta Rn
#> D_P  0.006489  1.425087  0.154132 -0.975352  1
#> 
#> Joint Wald statistic 2.031 , p-value 0.1541

qtest(univariate$Ret, univariate$D_P)
#>  CI for beta: [ -0.001058309 , -0.004163946 ]
#>  CI for rho: [ 0.9997793 , 1.003601 ]
#>  delta =  -0.9771747
```

And the multivariate estimation, for one or multiple horizons:

``` r
ivx(Ret ~ D_P + TBL, data = multivariate) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ D_P + TBL, data = multivariate, horizon = 1)
#> 
#> Coefficients:
#>      Estimate  Wald Ind Pr(> chi)     delta    Rn
#> D_P  0.006145  1.348538  0.177485 -0.975575 1.000
#> TBL -0.080717 -1.398871  0.161852 -0.061043 0.997
#> 
#> Joint Wald statistic 3.644 , p-value 0.1617

ivx(Ret ~ D_P + TBL, data = multivariate, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ D_P + TBL, data = multivariate, horizon = 4)
#> 
#> Coefficients:
#>      Estimate  Wald Ind Pr(> chi)     delta    Rn
#> D_P  0.006579  1.429957  0.152729 -0.975575 1.000
#> TBL -0.073549 -1.262909  0.206622 -0.061043 0.997
#> 
#> Joint Wald statistic 3.527 , p-value 0.1715
```


# ivx: Robust Econometric Inference

[![CRAN
status](https://www.r-pkg.org/badges/version/ivx)](https://cran.r-project.org/package=ivx)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/kvasilopoulos/ivx.svg?branch=master)](https://travis-ci.org/kvasilopoulos/ivx)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/kvasilopoulos/ivx?branch=master&svg=true)](https://ci.appveyor.com/project/kvasilopoulos/ivx)

Conducting inference on the regression coefficient with a highly
persistent regressor, using the the IVX approach for univariate
(Magdalinos and Phillips, 2009b) and multivariate (Kostakis Magdalinos
and Stamatogiannis, 2015) regression. Testing can be used for
long-horizon predictability as well.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kvasilopoulos/ivx")
```

## Usage

Load ivx and dplyr for data managment

``` r
library(ivx)
```

This is a basic example, lets load the data first:

``` r
# Monthly data from Kostakis et al (2014)
monthly %>%
  names()
#>  [1] "Date" "DE"   "LTY"  "DY"   "DP"   "TBL"  "EP"   "BM"   "INF"  "DFY" 
#> [11] "NTIS" "TMS"  "Ret"
```

## Univariate

And then do the univariate estimation:

``` r
ivx(Ret ~ DP, data = monthly) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = monthly, horizon = 1)
#> 
#> Coefficients:
#>    Estimate Wald Ind Pr(> chi)
#> DP 0.006489    2.018     0.155
#> 
#> Joint Wald statistic:  2.018 on 1 DF, p-value 0.1554

ivx(Ret ~ DP, data = monthly, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = monthly, horizon = 4)
#> 
#> Coefficients:
#>    Estimate Wald Ind Pr(> chi)
#> DP 0.006931    2.257     0.133
#> 
#> Joint Wald statistic:  2.257 on 1 DF, p-value 0.133
```

## Multivariate

And the multivariate estimation, for one or multiple horizons:

``` r
ivx(Ret ~ DP + TBL, data = monthly) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = monthly, horizon = 1)
#> 
#> Coefficients:
#>      Estimate Wald Ind Pr(> chi)
#> DP   0.006145    1.809     0.179
#> TBL -0.080717    1.956     0.162
#> 
#> Joint Wald statistic:  3.628 on 2 DF, p-value 0.163

ivx(Ret ~ DP + TBL, data = monthly, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = monthly, horizon = 4)
#> 
#> Coefficients:
#>      Estimate Wald Ind Pr(> chi)
#> DP   0.006579    2.034     0.154
#> TBL -0.073549    1.594     0.207
#> 
#> Joint Wald statistic:   3.51 on 2 DF, p-value 0.1729
```

-----

Please note that the ‘ivx’ project is released with a [Contributor Code
of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

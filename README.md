
# ivx: Robust Econometric Inference <img src='man/figures/logo.png' align="right" height="136.5" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ivx)](https://cran.r-project.org/package=ivx)
[![DOI](https://zenodo.org/badge/137785074.svg)](https://zenodo.org/badge/latestdoi/137785074)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R build
status](https://github.com/kvasilopoulos/ivx/workflows/R-CMD-check/badge.svg)](https://github.com/kvasilopoulos/ivx/actions)
[![codecov](https://codecov.io/gh/kvasilopoulos/ivx/branch/master/graph/badge.svg)](https://codecov.io/gh/kvasilopoulos/ivx)
<!-- badges: end -->

Drawing statistical inference on the coefficients of a short- or
long-horizon predictive regression with persistent regressors by using
the IVX method of [Magdalinos and Phillips
(2009)](https://doi.org/10.1017/S0266466608090154) and [Kostakis,
Magdalinos and Stamatogiannis
(2015)](https://doi.org/10.1093/rfs/hhu139).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# Install release version from CRAN
install.packages("ivx")


# install.packages("devtools")
devtools::install_github("kvasilopoulos/ivx")
```

## Usage

``` r
library(ivx)
library(magrittr)
```

This is a basic example, lets load the data first:

``` r
# Monthly data from Kostakis et al (2014)
kms %>%
  names()
#>  [1] "Date" "DE"   "LTY"  "DY"   "DP"   "TBL"  "EP"   "BM"   "INF"  "DFY"  "NTIS" "TMS" 
#> [13] "Ret"
```

## Univariate

And then do the univariate estimation:

``` r
ivx(Ret ~ DP, data = kms) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, horizon = 1)
#> 
#> Coefficients:
#>    Estimate Wald Ind Pr(> chi)
#> DP  0.00649     2.03      0.15
#> 
#> Joint Wald statistic:  2.03 on 1 DF, p-value 0.154
#> Multiple R-squared:  0.00284,    Adjusted R-squared:  0.00188

ivx(Ret ~ DP, data = kms, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, horizon = 4)
#> 
#> Coefficients:
#>    Estimate Wald Ind Pr(> chi)
#> DP  0.00693     2.27      0.13
#> 
#> Joint Wald statistic:  2.27 on 1 DF, p-value 0.132
#> Multiple R-squared:  0.0117, Adjusted R-squared:  0.0136
```

## Multivariate

And the multivariate estimation, for one or multiple horizons:

``` r
ivx(Ret ~ DP + TBL, data = kms) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 1)
#> 
#> Coefficients:
#>     Estimate Wald Ind Pr(> chi)
#> DP   0.00615     1.82      0.18
#> TBL -0.08072     1.96      0.16
#> 
#> Joint Wald statistic:  3.64 on 2 DF, p-value 0.162
#> Multiple R-squared:  0.00497,    Adjusted R-squared:  0.00304

ivx(Ret ~ DP + TBL, data = kms, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 4)
#> 
#> Coefficients:
#>     Estimate Wald Ind Pr(> chi)
#> DP   0.00658     2.04      0.15
#> TBL -0.07355     1.59      0.21
#> 
#> Joint Wald statistic:  3.53 on 2 DF, p-value 0.171
#> Multiple R-squared:  0.018,  Adjusted R-squared:  0.0189
```

## Yang et al. (2020) IVX-AR methodology

``` r
ivx_ar(hpi ~ cpi, data = ylpc) %>% 
  summary()
#> 
#> Call:
#> ivx_ar(formula = hpi ~ cpi, data = ylpc, horizon = 1)
#> 
#> Auto (bic) with AR terms q = 4
#> 
#> Coefficients:
#>      Estimate Wald Ind Pr(> chi)  
#> cpi -0.000177     4.33     0.038 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Joint Wald statistic:  4.33 on 1 DF, p-value 0.0375
#> Multiple R-squared:  0.0272, Adjusted R-squared:  0.0214
#> Wald AR statistic:  132 on 4 DF, p-value <0.0000000000000002
```

<!--
#### To-do  

* The Bonferroni method 
  - Cavanagh et al (1995)   
  - Campbell and Yogo (2006)    
* A conditional likelihood approach 
  - Jansson and Moreira (2006)  
* A control function approach   
  - Elliot (2001)
-->

-----

Please note that the ‘ivx’ project is released with a [Contributor Code
of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

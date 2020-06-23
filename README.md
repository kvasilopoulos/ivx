
# ivx: Robust Econometric Inference <img src='man/figures/logo.png' align="right" height="136.5" />

[![CRAN
status](https://www.r-pkg.org/badges/version/ivx)](https://cran.r-project.org/package=ivx)
[![DOI](https://zenodo.org/badge/137785074.svg)](https://zenodo.org/badge/latestdoi/137785074)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/kvasilopoulos/ivx.svg?branch=master)](https://travis-ci.org/kvasilopoulos/ivx)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/kvasilopoulos/ivx?branch=master&svg=true)](https://ci.appveyor.com/project/kvasilopoulos/ivx)
[![codecov](https://codecov.io/gh/kvasilopoulos/ivx/branch/master/graph/badge.svg)](https://codecov.io/gh/kvasilopoulos/ivx)

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
#>  [1] "Date" "DE"   "LTY"  "DY"   "DP"   "TBL"  "EP"   "BM"   "INF"  "DFY" 
#> [11] "NTIS" "TMS"  "Ret"
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
#> DP 0.006489    2.031     0.154
#> 
#> Joint Wald statistic:  2.031 on 1 DF, p-value 0.1541
#> Multiple R-squared:  0.002844,   Adjusted R-squared:  0.001877

ivx(Ret ~ DP, data = kms, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, horizon = 4)
#> 
#> Coefficients:
#>    Estimate Wald Ind Pr(> chi)
#> DP 0.006931    2.271     0.132
#> 
#> Joint Wald statistic:  2.271 on 1 DF, p-value 0.1318
#> Multiple R-squared:  0.01167,    Adjusted R-squared:  0.01358
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
#>      Estimate Wald Ind Pr(> chi)
#> DP   0.006145    1.819     0.177
#> TBL -0.080717    1.957     0.162
#> 
#> Joint Wald statistic:  3.644 on 2 DF, p-value 0.1617
#> Multiple R-squared:  0.004968,   Adjusted R-squared:  0.003036

ivx(Ret ~ DP + TBL, data = kms, horizon = 4) %>% 
  summary()
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 4)
#> 
#> Coefficients:
#>      Estimate Wald Ind Pr(> chi)
#> DP   0.006579    2.045     0.153
#> TBL -0.073549    1.595     0.207
#> 
#> Joint Wald statistic:  3.527 on 2 DF, p-value 0.1715
#> Multiple R-squared:  0.018,  Adjusted R-squared:  0.01895
```

## Yang et al. (2020) IVX-AR methodology

``` r
ivx_ar(hpi ~ cpi, data = ylpc)
#> 
#> Call:
#> ivx_ar(formula = hpi ~ cpi, data = ylpc, horizon = 1)
#> 
#> Auto (bic) with AR terms q = 4
#> 
#> Coefficients:
#>        cpi  
#> -0.0001775
```

#### To-do

  - The Bonferroni method
      - Cavanagh et al (1995)  
      - Campbell and Yogo (2006)  
  - A conditional likelihood approach
      - Jansson and Moreira (2006)  
  - A control function approach
      - Elliot (2001)

-----

Please note that the ‘ivx’ project is released with a [Contributor Code
of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

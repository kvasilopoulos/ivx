
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ivx: Predictive Regressions

[![CRAN
status](https://www.r-pkg.org/badges/version/ivx)](https://cran.r-project.org/package=ivx)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/kvasilopoulos/ivx.svg?branch=master)](https://travis-ci.org/kvasilopoulos/ivx)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/kvasilopoulos/ivx?branch=master&svg=true)](https://ci.appveyor.com/project/kvasilopoulos/ivx)

The goal of ivx is to offer robust econometric inference for predictive
regressions.

![first
equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20y_t%20%26%3D%20%5Cbeta%20y_%7Bt-1%7D%20+%20%5Cepsilon_t%5C%5C%20x_t%20%26%3D%20%5Crho%20x_%7Bt-1%7D%20+%20u_t%20%5Cend%7Balign*%7D)

### Solution methods

#### The workforce of the model hence the name is:

  - The IVX approach
      - Magdalinos and Phillips (2009b)
      - Kostakis Magdalinos and Stamatogiannis (2015)

#### The rest

  - The Bonferroni method
      - Cavanagh et al (1995)
      - Campbell and Yogo (2006)
  - A conditional likelihood approach
      - Jansson and Moreira (2006)
  - A control function approach
      - Elliot (2001)

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
#>     Estimate Wald Ind Pr(> chi)
#> D_P 0.006489    2.018     0.155
#> 
#> Joint Wald statistic:  2.018 on 1 DF, p-value 0.1554

qtest(univariate$Ret, univariate$D_P)
#>  CI for beta: [ -0.001058309 , -0.004163946 ]
#>  CI for rho: [ 0.9997793 , 1.003601 ]
#>  delta =  -0.9771747
```

And the multivariate estimation, for one or multiple horizons:

Please note that the ‘ivx’ project is released with a [Contributor Code
of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

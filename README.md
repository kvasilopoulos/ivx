
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
#> Warning in cor(cbind(epshat, u)): the standard deviation is zero
str(est)
#> List of 9
#>  $ Aivx   : num [1, 1:2] 0 -0.00329
#>  $ Wivx   : num [1:2, 1] 0.393 0.822
#>  $ WivxInd: num [1:2, 1:2] NaN NaN -0.627 0.531
#>  $ Q      : num [1:2, 1:2] 0.00 0.00 2.75e-05 2.75e-05
#>  $        : num [1:3, 1:3] 1 NA -0.0671 NA 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "" "(Intercept)" "X"
#>   .. ..$ : chr [1:3] "" "(Intercept)" "X"
#>  $ xlevels: Named list()
#>  $ call   : language ivx(formula = Y ~ X, data = reg)
#>  $ terms  :Classes 'terms', 'formula'  language Y ~ X
#>   .. ..- attr(*, "variables")= language list(Y, X)
#>   .. ..- attr(*, "factors")= int [1:2, 1] 0 1
#>   .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. ..$ : chr [1:2] "Y" "X"
#>   .. .. .. ..$ : chr "X"
#>   .. ..- attr(*, "term.labels")= chr "X"
#>   .. ..- attr(*, "order")= int 1
#>   .. ..- attr(*, "intercept")= int 1
#>   .. ..- attr(*, "response")= int 1
#>   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>   .. ..- attr(*, "predvars")= language list(Y, X)
#>   .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
#>   .. .. ..- attr(*, "names")= chr [1:2] "Y" "X"
#>  $ model  :'data.frame': 1033 obs. of  2 variables:
#>   ..$ Y: num [1:1033] 0.02296 -0.00522 0.04204 0.00457 0.01052 ...
#>   ..$ X: num [1:1033] -0.586 -0.568 -0.549 -0.531 -0.513 ...
#>   ..- attr(*, "terms")=Classes 'terms', 'formula'  language Y ~ X
#>   .. .. ..- attr(*, "variables")= language list(Y, X)
#>   .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
#>   .. .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. .. ..$ : chr [1:2] "Y" "X"
#>   .. .. .. .. ..$ : chr "X"
#>   .. .. ..- attr(*, "term.labels")= chr "X"
#>   .. .. ..- attr(*, "order")= int 1
#>   .. .. ..- attr(*, "intercept")= int 1
#>   .. .. ..- attr(*, "response")= int 1
#>   .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>   .. .. ..- attr(*, "predvars")= language list(Y, X)
#>   .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
#>   .. .. .. ..- attr(*, "names")= chr [1:2] "Y" "X"
#>  - attr(*, "class")= chr "ivx"
```

# Autocorrelation tests

Autocorrelation tests

## Usage

``` r
ac_test(x, lag_max = 5)
```

## Arguments

- x:

  the residuals or an `ivx` object.

- lag_max:

  the maximum length of lags.

## Examples

``` r
obj <- ivx(hpi ~ cpi + def + int + log(res), data = ylpc)
lmtest::bgtest(hpi ~ cpi + def + int + log(res), data = ylpc)
#> Error in loadNamespace(x): there is no package called 'lmtest'
ac_test(obj, 5)
#>  Lag     Wald LjungBox BoxPierce BreuschGodfrey
#>    1 27.68*** 61.86***   60.8***       61.87***
#>    2 48.65*** 91.56***  89.82***       62.89***
#>    3  47.7*** 145.2***  141.9***       84.85***
#>    4 126.9*** 206.8***  201.4***       90.76***
#>    5 125.4*** 236.4***  229.9***       90.81***
```

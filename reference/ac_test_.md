# Tests for autocorrelation

- `ac_test`: For all tests

- `ac_test_wald`: Wald test

- `ac_test_lb`: Ljung-Box

- `ac_test_bp`: Box-Pierce

- `ac_test_bg`: Breusch-Godfrey

## Usage

``` r
ac_test_wald(x, lag)

ac_test_lb(x, lag)

ac_test_bp(x, lag)

ac_test_bg(x, order, type, fill)
```

## Arguments

- x:

  an `ivx` model or a `numeric vector`, usually the residuals from an
  ols regression.

- lag:

  the number of lags.

- order:

  lag TODO

- type:

  the type of test statistic to be returned. Either "Chisq" for the
  Chi-squared test statistic or "F" for the F test statistic.

- fill:

  starting values for the lagged residuals in the auxiliary regression.
  By default 0 but can also be set to NA.

## Value

a numeric scalar or numeric vector.

## Details

If p-value \< 0.051: You can reject the null hypothesis assuming a 5%
chance of making a mistake. So you can assume that your values are
showing dependence on each other.

## See also

`Box.test`
[`lmtest::bgtest`](https://rdrr.io/pkg/lmtest/man/bgtest.html)

## Examples

``` r

mdl <- ivx(hpi ~ cpi + inv, data = ylpc)
ac_test_wald(mdl)
#>  Lag     Wald
#>    1 17.39***

ac_test(mdl)
#>  Lag     Wald LjungBox BoxPierce BreuschGodfrey
#>    1 17.39*** 42.56***  41.83***       46.73***
#>    2 45.31*** 57.45***  56.38***       46.94***
#>    3 47.96*** 99.86***  97.57***       70.99***
#>    4 82.82*** 155.5***  151.3***       82.14***
#>    5 91.79*** 176.6***  171.6***       82.31***
```

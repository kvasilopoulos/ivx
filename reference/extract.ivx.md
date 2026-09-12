# `extract` method for `ivx` objects

`extract` method for `ivx` objects

## Usage

``` r
extract.ivx(
  model,
  include.wald = TRUE,
  include.nobs = TRUE,
  include.aic = FALSE,
  include.bic = FALSE,
  include.rsquared = FALSE,
  include.adjrs = FALSE,
  ...
)

extract.ivx_ar(
  model,
  include.wald = TRUE,
  include.nobs = TRUE,
  include.aic = FALSE,
  include.bic = FALSE,
  include.rsquared = FALSE,
  include.adjrs = FALSE,
  ...
)
```

## Arguments

- model:

  A statistical model object.

- include.wald:

  Report the Wald statistic.

- include.nobs:

  Report the number of observations in the GOF block?

- include.aic:

  Report Akaike's Information Criterion (AIC) in the GOF block?

- include.bic:

  Report the Bayesian Information Criterion (BIC) in the GOF block?

- include.rsquared:

  Report the R-squared.

- include.adjrs:

  Report the Adjusted R-squared.

- ...:

  Custom parameters, which are handed over to subroutines. Currently not
  in use.

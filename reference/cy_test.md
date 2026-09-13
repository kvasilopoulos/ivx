# Bonferroni Q-Test of Campbell and Yogo (2006)

`cy_test` implements the Bonferroni Q-test of Campbell and Yogo (2006),
the standard feasible version of the sup-bound / Bonferroni approach of
Cavanagh, Elliott and Stock (1995). For a single predictor \\x_t =
\gamma + \rho x\_{t-1} + v_t\\ (AR(p) dynamics allowed, Appendix A), the
procedure is:

1.  estimate the innovation correlation \\\delta\\ between the
    predictive-regression residual and the ADF innovation of the
    predictor (lag length by BIC, \\p \in \[1, p\_{max}\]\\);

2.  compute the DF-GLS statistic of Elliott, Rothenberg and Stock (1996)
    and invert its local-to-unity null distribution (Stock, 1991) into a
    confidence interval \\\[\underline c, \bar c\]\\ for \\c = T(\rho -
    1)\\, with lower and upper levels \\\underline\alpha_1(\delta)\\,
    \\\bar\alpha_1(\delta)\\ from the paper's Table 2, which refine the
    Bonferroni bound so that the one-sided test has size 5\\

3.  for \\\rho\\ at each end of the interval compute the Q-estimate
    \\\hat\beta(\rho)\\ of eq. (25) - OLS of \\y_t\\ on the demeaned
    \\x\_{t-1}\\ after subtracting \\(\sigma\_{ue}/\sigma_e\omega)(x_t -
    \rho x\_{t-1})\\, with the Phillips-Perron-type correction
    \\\tfrac{T}{2}(\sigma\_{ue}/\sigma_e\omega)(\omega^2 - \sigma_v^2)\\
    when \\p \> 1\\ - and its standard error \\\sigma_u
    (1-\delta^2)^{1/2}/(\sum x^{\mu 2}\_{t-1})^{1/2}\\;

4.  the 90\\ \\\[\hat\beta(\bar\rho) - 1.645\\se,\\
    \hat\beta(\underline\rho) + 1.645\\se\]\\ (eq. 17); the 5\\ bound is
    positive and \\\beta \ge 0\\ if the upper bound is negative.

The DF-GLS null quantiles are tabulated by simulation for \\c \in
\[-100, 10\]\\ (see `data-raw/dfgls-quantiles.R`); a statistic outside
the tabulated range is clamped to the boundary, which for very negative
values (a clearly stationary predictor) makes the interval for \\c\\
start at \\-100\\. Table 2 is given for \\\delta \< 0\\; for
\\\hat\delta \> 0\\ the predictor is sign-flipped, which flips \\\beta\\
and the alternative, and the results are mapped back.

## Usage

``` r
cy_test(formula, data, lag_max = NULL, na.action)

cy_test_fit(y, x, lag_max = NULL)

# S3 method for class 'cy_test'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted.

- data:

  n optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in data,
  the variables are taken from environment(formula), typically the
  environment from which lm is called.

- lag_max:

  maximum ADF lag order for the BIC search; the default is \\\lfloor 12
  (T/100)^{1/4} \rfloor\\ lagged differences.

- na.action:

  a function which indicates what should happen when the data contain
  NAs. The default is set by the na.action setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

- y:

  response vector.

- x:

  an object of class "cy_test".

- digits:

  minimal number of significant digits.

- ...:

  unused.

## Value

an object of class "cy_test": a list with `ci` (the 90\\ interval for
\\\beta\\), `reject` (named logical: `greater`, `less`), `estimate` (OLS
slope), `beta_rho` (\\\hat\beta\\ at \\\underline\rho\\ and
\\\bar\rho\\), `se`, `delta`, `dfgls`, `c_ci`, `rho_ci`, `alpha1`, `lag`
and `n`.

## References

Campbell, J. Y., & Yogo, M. (2006). Efficient tests of stock return
predictability. Journal of Financial Economics, 81(1), 27-60.

Cavanagh, C. L., Elliott, G., & Stock, J. H. (1995). Inference in models
with nearly integrated regressors. Econometric Theory, 11(5), 1131-1147.

Elliott, G., Rothenberg, T. J., & Stock, J. H. (1996). Efficient tests
for an autoregressive unit root. Econometrica, 64(4), 813-836.

Stock, J. H. (1991). Confidence intervals for the largest autoregressive
root in U.S. macroeconomic time series. Journal of Monetary Economics,
28(3), 435-459.

## Examples

``` r
cy_test(Ret ~ DP, data = kms)
#> 
#> Call:
#> cy_test(formula = Ret ~ DP, data = kms)
#> 
#> Bonferroni Q-test (Campbell & Yogo, 2006)
#> 
#> delta = -0.972, DF-GLS = -1.468 (p = 2), CI for c at levels (0.055, 0.082): [-9.319, 1.045], rho: [0.991, 1.001]
#> OLS slope = 0.006193; Q-estimates at the ends of the rho interval: 0.009076, 0.0004903
#> 90% Bonferroni confidence interval for beta: [-0.000973, 0.01054]
#> 5% one-sided Q-tests: H1 beta > 0 do not reject H0; H1 beta < 0 do not reject H0
#> 
```

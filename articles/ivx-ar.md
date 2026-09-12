# IVX-AR: serially correlated errors

``` r

library(ivx)
```

[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
implements the IVX-AR procedure of Yang, Long, Peng & Cai (2020), which
extends the Kostakis et al. (2015) test to predictive regressions whose
errors are serially correlated — the situation the authors document for
housing price index returns.

## Model

\\ y_t = \mu + \beta' x\_{t-1} + u_t, \qquad u_t = \sum\_{j=1}^{q}
\phi_j u\_{t-j} + \varepsilon_t, \\

with persistent predictors as in
[`vignette("ivx")`](https://kvasilopoulos.github.io/ivx/articles/ivx.md).
Serial correlation in \\u_t\\ invalidates the IVX Wald statistic’s
\\\chi^2\\ limit (and the Newey–West correction inside it is designed
for the predictor’s innovations, not for \\u_t\\). Applying the AR
filter \\\phi(L) = 1 - \sum_j \phi_j L^j\\ to both sides,

\\ \phi(L)\\ y_t = \phi(1)\mu + \beta' \phi(L)\\ x\_{t-1} +
\varepsilon_t, \\

gives a predictive regression with white-noise errors and the same slope
\\\beta\\, to which IVX applies. Since \\\phi\\ is unknown and its
estimate from the residuals is biased under persistence, the paper
profiles over it.

## Procedure

1.  Fit the plain IVX regression and fit an AR(\\q\\) model to its OLS
    residuals, choosing \\q\\ by an information criterion
    (`ar = "auto"`, `ar_ic`, `ar_max`) or fixing it (`ar = q`).
2.  For every point \\\phi\\ of a grid around the AR estimates
    (`ar_grid`, default \\\pm 0.3\\ in steps of 0.02 per coefficient),
    quasi-difference \\y_t\\ and \\x\_{t-1}\\ with \\\phi(L)\\ and refit
    IVX.
3.  Keep the grid point with the smallest residual variance; report that
    fit’s IVX Wald statistics.
4.  Additionally test \\H_0: \phi_1 = \dots = \phi_q = 0\\ with the Wald
    statistic `Wald_AR` (see
    [`ac_test_wald()`](https://kvasilopoulos.github.io/ivx/reference/ac_test_.md));
    `ar = 0` reduces to
    [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md).

``` r

m <- ivx_ar(hpi ~ log(res) + cpi, data = ylpc)
m
#> 
#> Call:
#> ivx_ar(formula = hpi ~ log(res) + cpi, data = ylpc, horizon = 1)
#> 
#> Lag Selection:
#> Auto (bic) with AR terms q = 4
#> 
#> Coefficients:
#>   log(res)         cpi  
#>  0.0018317  -0.0002078
summary(m)
#> 
#> Call:
#> ivx_ar(formula = hpi ~ log(res) + cpi, data = ylpc, horizon = 1)
#> 
#> Auto () with AR terms q = 4
#> 
#> Coefficients:
#>            Estimate Std. Error t value Wald Ind Pr(> chi)  
#> log(res)  0.0018317  0.0012533   1.461    2.136    0.1439  
#> cpi      -0.0002078  0.0001075  -1.932    3.732    0.0534 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Joint Wald statistic:  4.316 on 2 DF, p-value 0.1156
#> Multiple R-squared:  0.03437,    Adjusted R-squared:  0.02281
#> Wald AR statistic: 133.6 on 4 DF, p-value < 2.2e-16
```

The estimated AR coefficients of the selected grid point are in
`m$coefficients_ar`, and the automatic order selection can be replaced:

``` r

m$coefficients_ar
#>        ar1        ar2        ar3        ar4 
#>  0.4062933 -0.1139760  0.3459206  0.2142420
coef(ivx_ar(hpi ~ log(res) + cpi, data = ylpc, ar = 1))
#>      log(res)           cpi 
#>  0.0010780144 -0.0001155392
```

The `ylpc` dataset is the authors’ quarterly US housing data.

## Caveats

- The grid search is over \\q\\ coefficients with the default 31 points
  each, so cost grows as \\31^q\\; keep `ar_max` small or fix `ar`.
- `ar = "forecast"` delegates order selection to
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
  and needs that package.
- The tuning and standard-error options of
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)
  (`beta`, `cz`, `bandwidth`, `robust`) are passed through to every grid
  refit.
- [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  does not accept `ivx_ar` objects: the bootstrap schemes are defined
  for the untransformed regression.
- When the residuals show no serial correlation
  ([`ac_test()`](https://kvasilopoulos.github.io/ivx/reference/ac_test.md)
  on an `ivx` fit is a quick check)
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) is the
  better choice; the profiling adds noise.

## References

- Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the
  predictability of US housing price index returns based on an IVX-AR
  model. *Journal of the American Statistical Association*, 115(532),
  1598–1619.
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). Robust
  econometric inference for stock return predictability. *Review of
  Financial Studies*, 28(5), 1506–1553.

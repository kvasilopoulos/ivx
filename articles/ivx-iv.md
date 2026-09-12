# IV tests with combined instruments

``` r

library(ivx)
```

[`ivx_iv()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md)
implements the instrumental-variable predictability tests of Breitung &
Demetrescu (2015). IVX is one member of the family: any instrument
\\z\_{t-1}\\ that is less persistent than \\x\_{t-1}\\ but correlated
with it gives a 2SLS \\t\\-ratio with a standard normal null limit
whatever the local-to-unity parameter \\c\\.

## Instruments

Type-I instruments are functions of the predictor:

- `"frac"`: the truncated fractional difference \\\Delta\_+^{d}
  x\_{t-1}\\, \\d = 1/2\\ by default — \\I(1/2)\\ when \\x\\ is
  \\I(1)\\, still informative when \\x\\ is stationary;
- `"diff"`: the long difference \\x\_{t-1} - x\_{t-1-k_T}\\ with \\k_T =
  \lfloor 0.2\\T^{0.85} \rfloor\\;
- the IVX filter of Kostakis et al. (2015) is the paper’s “mild
  integration” case
  ([`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)).

Type-II instruments are deterministic and only correlate with a
near-integrated predictor: `"sin"` uses \\\sin(\pi t/T)\\ (with \\K\\
predictors, frequencies \\k = 1, \dots, K\\). They are weak when \\x\\
is stationary, so the paper recommends the over-identified combination
`"comb"` = sine + fractional: the 2SLS statistic is asymptotically
driven by whichever instrument is informative (Theorem 3), so it keeps
power in both regimes.

## The statistic

With \\Z\\ the instrument matrix (intercept included in both stages) and
\\\hat u_t\\ the OLS residuals of the predictive regression, the squared
\\t\\-ratio of eq. (12) is \\ t^2 = \frac{\big(x'Z(Z'Z)^{-1}Z'y\big)^2}
{x'Z(Z'Z)^{-1}\big(\sum_t z_t z_t' \hat u_t^2\big)(Z'Z)^{-1}Z'x} \\\to\\
\chi^2(1), \\ and its \\K\\-predictor Wald analogue is \\\chi^2(K)\\
(Theorem 7).

``` r

summary(ivx_iv(Ret ~ DP + TBL, data = kms))
#> 
#> Call:
#> ivx_iv(formula = Ret ~ DP + TBL, data = kms)
#> 
#> Coefficients:
#>       Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP   0.0006284  0.0061730   0.102    0.010     0.919
#> TBL -0.0626958  0.0822806  -0.762    0.581     0.446
#> (Eicker-White standard errors)
#> 
#> Joint Wald statistic:  0.9973 on 2 DF, p-value 0.6073
#> Multiple R-squared:  0.007404,   Adjusted R-squared:  0.00451
sapply(c("comb", "sin", "frac", "diff"),
       function(i) ivx_iv(Ret ~ DP, data = kms, instruments = i)$tstat)
#>    comb.DP     sin.DP    frac.DP    diff.DP 
#>  0.4900820  0.4655899 -0.7089655 -1.6445824
```

## Replication of Table 1

The paper’s design: \\T = 250\\, \\\mathrm{corr}(u, v) = 0.9\\, \\\beta
= b/T\\, two-sided tests at 10%. Rejection rates in 1000 replications
(paper values in brackets, 10 000 replications):

| \\\rho\\ | \\b\\ | comb        | sin         | frac        | diff        |
|----------|-------|-------------|-------------|-------------|-------------|
| 1.00     | 0     | 12.7 (11.2) | 11.5 (9.9)  | 12.2 (11.1) | 11.2 (12.5) |
| 1.00     | 10    | 52.6 (65.7) | 38.5 (61.4) | 41.5 (33.4) | 30.4 (33.6) |
| 1.00     | 20    | 88.5 (91.3) | 59.6 (79.5) | 81.0 (66.9) | 68.9 (66.1) |
| 0.96     | 0     | 11.6 (10.4) | 9.7 (9.9)   | 11.4 (9.4)  | 10.7 (10.5) |
| 0.96     | 10    | 42.6 (47.3) | 29.4 (34.4) | 31.4 (30.1) | 33.2 (33.4) |

Sizes and the ranking (combination best, sine strong only near the unit
root) reproduce; the split of power between the sine and fractional
instruments differs from the paper, which does not report the exact
demeaning and initialisation used for the deterministic instrument.

## Caveats

- Short horizon only, no bootstrap.
- The sign (Cauchy) instrument of the paper needs recursive demeaning of
  the regressor and forward demeaning of the response and is not
  included.
- Several type-II instruments per predictor are not allowed (Assumption
  5); the implementation uses one sine frequency per predictor.

## References

- Breitung, J., & Demetrescu, M. (2015). Instrumental variable and
  variable addition based inference in predictive regressions. *Journal
  of Econometrics*, 187(1), 358–375.

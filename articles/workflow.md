# Choosing a test: a decision workflow

``` r

library(ivx)
```

The package holds a dozen estimators and tests. Each vignette documents
one of them; this one is the map. It runs a predictive regression from
the plain
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) fit
through the diagnostics that decide which extension applies, and follows
every branch with the `kms` (monthly S&P 500 excess returns, Kostakis et
al. 2015) and `ylpc` (quarterly US housing returns, Yang et al. 2020)
data.

## The map

    ivx(y ~ x, data)                       KMS Wald test, the baseline
     │
     ├─ residuals serially correlated?  ac_test()          → ivx_ar()
     ├─ residuals heteroskedastic?      ac_test(res^2)     → ivx(robust = TRUE)
     ├─ |delta| near 1 and x near a unit root, short sample, one-sided test?
     │                                   delta(), hlt_test  → ivx_boot(), ivx_ra()
     ├─ weak signal, power a concern?                      → ivx(lag_y = TRUE)
     ├─ K-period cumulative return?                        → ivx(horizon = K), ivx_ra(horizon = K)
     ├─ several responses, cross-equation restrictions?    → ivx_sys()
     ├─ tail or median predictability?                     → ivx_qr(), ivx_qr_boot()
     ├─ predictability confined to episodes?               → ivx_episodic(), rolling ivx()
     └─ single predictor, want a non-IVX cross-check?      → hlt_test(), cy_test(), arm(), el_test(), ivx_iv()

The branches are not exclusive: the worked example below ends up
combining three of them.

## Step 1: baseline fit and three diagnostics

``` r

m <- ivx(Ret ~ DP + TBL, data = kms)
summary(m)
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 1)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP   0.006145   0.004557   1.349    1.819     0.177
#> TBL -0.080717   0.057701  -1.399    1.957     0.162
#> 
#> Joint Wald statistic:  3.644 on 2 DF, p-value 0.1617
#> Multiple R-squared:  0.004968,   Adjusted R-squared:  0.003036
```

The IVX Wald test is valid for any degree of persistence in \\x_t\\, so
the first question is not “is \\x_t\\ persistent?” but “is anything else
wrong with the errors, and is the sample short enough for finite-sample
distortions to matter?”. Three checks answer it.

**Persistence and endogeneity.** The lag-1 autocorrelation of each
predictor and the correlation \\\delta\\ between the return innovation
and the predictor innovation:

``` r

sapply(kms[c("DP", "TBL")], function(x) cor(head(x, -1), tail(x, -1)))
#>        DP       TBL 
#> 0.9923205 0.9926221
delta(m)
#> [1] -0.9755750 -0.0610432
```

Both predictors are near unit roots. `DP` is also strongly endogenous,
\\\delta \approx -0.98\\, the combination under which the asymptotic
test over-rejects in samples of a few hundred observations
(Hosseinkouchack & Demetrescu 2021; Demetrescu et al. 2023). `TBL` is
persistent but nearly exogenous, so its asymptotic p-value can be taken
at face value.
[`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md)
classifies persistence formally through the ADF normalised-bias
statistic (see its `adf` and `test` components) if a rule is preferred
to eyeballing.

**Serial correlation in the errors.**
[`ac_test()`](https://kvasilopoulos.github.io/ivx/reference/ac_test.md)
on the fit reports four tests at lags 1 to `lag_max`. The `Wald` column
is the regression-based test of Yang et al. (2020), the one
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
reports as `Wald_AR`; the Box tests are portmanteau tests on the
residual autocorrelations.

``` r

ac_test(m, 4)
#>  Lag    Wald LjungBox BoxPierce BreuschGodfrey
#>    1 0.8348   8.63***  8.605***       9.497***
#>    2  2.871   8.764**   8.739**       9.836***
#>    3  4.638  17.23***  17.17***       17.49***
#>    4  7.181  18.66***  18.59***       20.48***
```

The Box and Breusch–Godfrey statistics are significant, the Wald test is
not. The first-order autocorrelation of the residuals is 0.092: the mild
low-order dependence typical of monthly returns, which is not what
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md) is
for. Compare the housing data:

``` r

mh <- ivx(hpi ~ log(res) + cpi, data = ylpc)
ac_test(mh, 4)
#>  Lag     Wald LjungBox BoxPierce BreuschGodfrey
#>    1 31.15*** 65.89***  64.76***       65.28***
#>    2 38.57*** 100.4***  98.51***       66.64***
#>    3 66.62*** 160.5***  156.8***       88.84***
#>    4 133.6***   226***  220.1***       93.28***
```

Every statistic is an order of magnitude larger and the Wald test
rejects at lag 1. That is the
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
case (Step 2).

**Heteroskedasticity.** The same portmanteau tests on the squared
residuals are the McLeod–Li test for ARCH effects:

``` r

ac_test(residuals(m)^2, 4)
#>  Lag     Wald LjungBox BoxPierce
#>    1   1.407  58.23***  58.06***
#>    2 13.22*** 110.5***  110.1***
#>    3  23.1*** 205.5***  204.7***
#>    4 23.23***   240***  238.9***
```

Strong volatility clustering, as expected for returns. The Eicker–White
standard errors of Step 3 handle it.

## Step 2: serially correlated errors → `ivx_ar()`

Serial correlation in \\u_t\\ breaks the \\\chi^2\\ limit of the IVX
Wald statistic.
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
fits an AR(\\q\\) to the residuals, quasi-differences the regression
over a grid of AR coefficients and reports the IVX statistics of the fit
with the smallest residual variance
([`vignette("ivx-ar")`](https://kvasilopoulos.github.io/ivx/articles/ivx-ar.md)).

``` r

ma <- ivx_ar(hpi ~ log(res) + cpi, data = ylpc)
summary(ma)
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

`Wald_AR` is the test of \\\phi_1 = \dots = \phi_q = 0\\; when it does
not reject, the profiling adds noise and plain
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) is the
better choice.
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
passes `robust`, `beta`, `cz` and `bandwidth` through to every grid
refit, but
[`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
does not accept its output: the bootstrap is defined for the
untransformed regression.

## Step 3: heteroskedasticity → `robust = TRUE`

Demetrescu, Georgiev, Rodrigues & Taylor (2023) show the IVX statistics
keep their limits under conditional and unconditional heteroskedasticity
if the homoskedastic \\\hat\sigma^2_u \sum z z'\\ is replaced by the
Eicker–White \\\sum z\_{t-1} z\_{t-1}' \hat u_t^2\\. There is no cost
when the errors are homoskedastic, so for return data this is the
sensible default.

``` r

mr <- ivx(Ret ~ DP + TBL, data = kms, robust = TRUE)
coef(summary(mr))
#>         Estimate  Std. Error   t value Wald Ind Pr(> chi)
#> DP   0.006145163 0.004792361  1.282283 1.644250 0.1997434
#> TBL -0.080716672 0.057249850 -1.409902 1.987823 0.1585687
```

`robust = TRUE` is defined for `horizon = 1`;
[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
uses heteroskedasticity-robust standard errors at every horizon.

## Step 4: finite-sample size → `ivx_boot()` and `ivx_ra()`

For a strongly endogenous, near-unit-root predictor such as `DP`, the
asymptotic test is oversized in finite samples, more so one-sided. Two
remedies, in order of preference:

1.  **Bootstrap the statistic.** The residual wild bootstrap of
    DGRT (2023) refits the predictor’s autoregression, resamples \\(u_t,
    v_t)\\ with a common multiplier so \\\delta\\ is preserved, and
    recomputes the IVX statistics
    ([`vignette("robust-inference")`](https://kvasilopoulos.github.io/ivx/articles/robust-inference.md)).

``` r

mb <- ivx_boot(m, B = 499, type = "rwb", seed = 1)
mb
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 1)
#> 
#> Residual wild bootstrap, B = 499
#> 
#> Coefficients (bootstrap p-values):
#>      Estimate t value Wald Ind Pr(> chi) Pr(t < 0) Pr(t > 0)
#> DP   0.006145   1.349    1.819    0.3768    0.6693     0.331
#> TBL -0.080717  -1.399    1.957    0.2084    0.1002     0.900
#> 
#> Joint Wald statistic: 3.644, bootstrap p-value 0.3186
```

2.  **Bias-reduce the estimator.**
    [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
    (Demetrescu & Rodrigues 2022) augments the regression with the
    residuals of the predictors’ autoregression, which removes the
    endogeneity-driven bias before instrumenting. Its standard errors
    are robust by construction.

``` r

mra <- ivx_ra(Ret ~ DP + TBL, data = kms)
coef(summary(mra))
#>         Estimate  Std. Error    t value  Wald Ind Pr(> chi)
#> DP  -0.001934204 0.005010361 -0.3860409 0.1490275 0.6994664
#> TBL -0.053983764 0.055313355 -0.9759626 0.9525030 0.3290830
```

The four p-values for the two-sided test of `DP` line up as the theory
predicts, with the asymptotic homoskedastic one the most optimistic:

``` r

p_dp <- c(
  asymptotic = coef(summary(m))["DP", "Pr(> chi)"],
  eicker_white = coef(summary(mr))["DP", "Pr(> chi)"],
  rwb_bootstrap = mb$p.value$Wald_Ind[["DP"]],
  residual_augmented = coef(summary(mra))["DP", "Pr(> chi)"]
)
round(p_dp, 3)
#>         asymptotic       eicker_white      rwb_bootstrap residual_augmented 
#>              0.177              0.200              0.377              0.699
```

Use `cores > 1` for the bootstrap with \\B\\ in the thousands; use
`ivx_boot(type = "frwb")` if the predictor’s autoregression is unstable
(very short samples, structural breaks in the predictors).

## Step 5: power → `lag_y = TRUE`

The IVX instrument is less persistent than \\x_t\\, so the test gives up
local power against OLS exactly when \\x_t\\ is near-integrated and
endogenous. Demetrescu (2014) recovers part of it by adding \\y\_{t-1}\\
to the regression, instrumented by itself. This is a modelling choice
rather than a diagnostic outcome: it costs nothing asymptotically when
the extra term is irrelevant and can raise power substantially when the
signal is weak.

``` r

coef(summary(ivx(Ret ~ DP + TBL, data = kms, lag_y = TRUE)))
#>           Estimate  Std. Error   t value Wald Ind   Pr(> chi)
#> DP     0.007247938 0.004665519  1.553512 2.413398 0.120300968
#> TBL   -0.073288383 0.057667741 -1.270873 1.615119 0.203773782
#> y_lag  0.091950133 0.031084801  2.958042 8.750010 0.003096004
```

The joint Wald statistic still tests only the predictors.

## Step 6: what is the null?

The steps above fix the error structure. The remaining branches depend
on what is being tested.

**Long horizons.** For a \\K\\-period cumulative return, `horizon = K`
gives the Kostakis et al. (2023) statistic, which handles the overlap of
the sums through the instrument. `ivx_ra(horizon = K)` is the
transformed-regression alternative of Demetrescu, Rodrigues & Taylor
(2023), which avoids HAC estimation altogether. Report both when the
horizon is long relative to the sample: they use different information
and the overlap (\\K - 1\\ observations per statistic) is where
long-horizon tests misbehave.

``` r

c(kms = ivx(Ret ~ DP + TBL, data = kms, horizon = 12)$Wald_Joint,
  drt = ivx_ra(Ret ~ DP + TBL, data = kms, horizon = 12)$Wald_Joint)
#>      kms      drt 
#> 3.998219 1.200306
```

**Several responses.**
[`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
(Magdalinos 2022) estimates the system jointly and returns the full
covariance of \\\mathrm{vec}(A)\\, so restrictions across equations,
e.g. “does `DP` predict the return and the payout ratio equally?”, are
Wald tests on [`vcov()`](https://rdrr.io/r/stats/vcov.html)
([`vignette("ivx-sys")`](https://kvasilopoulos.github.io/ivx/articles/ivx-sys.md)).

**Quantiles.** If the question is about the tails (“does `DP` predict
large losses?”) rather than the mean,
[`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
(Lee 2016) runs the IVX test in a quantile regression;
[`ivx_qr_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_boot.md)
gives block-bootstrap intervals that are robust to the conditional
heteroskedasticity found in Step 1
([`vignette("ivx-qr")`](https://kvasilopoulos.github.io/ivx/articles/ivx-qr.md)).

``` r

q <- ivx_qr(Ret ~ DP + TBL, data = kms, tau = c(0.1, 0.5, 0.9))
sapply(q, function(fit) c(Wald_Joint = fit$Wald_Joint, p = 1 - pchisq(fit$Wald_Joint, 2)))
#>                   0.1        0.5          0.9
#> Wald_Joint 6.66613413 7.54143943 3.458798e+01
#> p          0.03568349 0.02303548 3.085432e-08
```

**Episodes.** A full-sample rejection can be driven by a short pocket of
predictability, and a full-sample non-rejection can hide one.
[`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md)
computes the IVX statistic over rolling, forward or backward windows
with the full-sample instrument and bootstraps the supremum of the
sequence
([`vignette("ivx-episodic")`](https://kvasilopoulos.github.io/ivx/articles/ivx-episodic.md));
a plain rolling
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) with
Bonferroni critical values is the bubble-detection variant of
[`vignette("rolling-ivx")`](https://kvasilopoulos.github.io/ivx/articles/rolling-ivx.md).

``` r

ivx_episodic(ivx(Ret ~ DP, data = kms), scheme = "rolling", window = 0.2, B = 199, seed = 1)
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, horizon = 1)
#> 
#> Subsample IVX tests, rolling scheme (window = 0.2), 827 windows
#> Fixed regressor wild bootstrap, B = 199
#> 
#>                         statistic bootstrap p
#> sup t   (H1: beta > 0)      2.921      0.1508
#> inf t   (H1: beta < 0)    -0.4166      1.0000
#> sup t^2 (H1: beta != 0)     8.533      0.2714
```

## Step 7: cross-check with a non-IVX test

With a single predictor, the benchmarks in
[`vignette("benchmarks")`](https://kvasilopoulos.github.io/ivx/articles/benchmarks.md)
answer the same question from a different direction and are cheap to
run. The hybrid test of Harvey, Leybourne & Taylor (2021) is the most
useful companion: it classifies the predictor’s persistence, picks the
test accordingly and reports the intermediate quantities.

``` r

hlt_test(Ret ~ DP, data = kms)
#> 
#> Call:
#> hlt_test(formula = Ret ~ DP, data = kms)
#> 
#> Hybrid predictability test (Harvey, Leybourne & Taylor, 2021)
#> 
#> Selected test: T~_con (quasi-GLS-demeaned t, conservative critical value)
#> statistic = 1.298, 5% critical value = 1.944 (alternative: beta > 0): do not reject
#> slope = 0.006172, t =  1.63, t (quasi-GLS) = 1.298, ADF = -5.97 (p = 21, cutoff -128.6), rho_xy = -0.9485
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
```

[`cy_test()`](https://kvasilopoulos.github.io/ivx/reference/cy_test.md)
(Campbell & Yogo 2006) gives a Bonferroni confidence interval for
\\\beta\\ that is exact for a near-unit-root predictor but conservative
otherwise;
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) and
[`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md)
are the reduced-bias OLS and empirical likelihood alternatives;
[`ivx_iv()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md)
(Breitung & Demetrescu 2015) replaces the IVX filter by sine or
fractional-difference instruments, a check that the result does not
hinge on the choice of \\\rho_n\\.

## The worked example, in one table

For `Ret ~ DP + TBL` the diagnostics said: near-unit-root predictors,
`DP` strongly endogenous, mild residual autocorrelation, strong ARCH.
The decisions that follow are Eicker–White standard errors, a bootstrap
p-value for the endogenous predictor, and
[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md) as
the bias-reduced cross-check; no
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md).
The p-values collected above, plus the joint test:

``` r

rbind(
  DP = round(p_dp, 3),
  joint = round(c(
    asymptotic = 1 - pchisq(m$Wald_Joint, 2),
    eicker_white = 1 - pchisq(mr$Wald_Joint, 2),
    rwb_bootstrap = mb$p.value$Wald_Joint,
    residual_augmented = 1 - pchisq(mra$Wald_Joint, 2)
  ), 3)
)
#>       asymptotic eicker_white rwb_bootstrap residual_augmented
#> DP         0.177        0.200         0.377              0.699
#> joint      0.162        0.235         0.319              0.491
```

Every robust version moves the p-value up, and the more the method does
about endogeneity the further it moves: the pattern DGRT (2023) document
for the dividend–price ratio. Here nothing rejects at any conventional
level, but with a borderline asymptotic p-value the ordering is what
decides the conclusion.

For `hpi ~ log(res) + cpi` the diagnostics said: strong residual
autocorrelation. The decision is
[`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md),
whose `Wald_AR` confirms the AR terms and whose IVX statistics replace
the plain ones.

## Summary

| Symptom | Diagnostic | Use |
|----|----|----|
| Serially correlated errors | `ac_test(fit)`, `Wald` column | [`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md) |
| Conditional heteroskedasticity | `ac_test(residuals(fit)^2)` | `ivx(robust = TRUE)` |
| \\\lvert\delta vert\\ near 1, \\x\\ near unit root, short sample, one-sided test | [`delta()`](https://kvasilopoulos.github.io/ivx/reference/delta.md), `hlt_test()$adf` | `ivx_boot(type = "rwb")`, [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md) |
| Weak signal | — | `ivx(lag_y = TRUE)` |
| Multi-period return | — | `ivx(horizon = K)`, `ivx_ra(horizon = K)` |
| Several responses | — | [`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md) |
| Tail predictability | — | [`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md), [`ivx_qr_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_boot.md) |
| Pockets of predictability | plot of rolling [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) t-ratios | [`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md) |
| Single predictor, second opinion | — | [`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md), [`cy_test()`](https://kvasilopoulos.github.io/ivx/reference/cy_test.md), [`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md), [`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md), [`ivx_iv()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md) |

## References

- Breitung, J., & Demetrescu, M. (2015). Instrumental variable and
  variable addition based inference in predictive regressions. *Journal
  of Econometrics*, 187(1), 358–375.
- Campbell, J. Y., & Yogo, M. (2006). Efficient tests of stock return
  predictability. *Journal of Financial Economics*, 81(1), 27–60.
- Demetrescu, M. (2014). Enhancing the local power of IVX-based tests in
  predictive regressions. *Economics Letters*, 124(2), 269–273.
- Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
  (2023). Extensions to IVX methods of inference for return
  predictability. *Journal of Econometrics*, 237(2), 105271.
- Demetrescu, M., & Rodrigues, P. M. M. (2022). Residual-augmented IVX
  predictive regression. *Journal of Econometrics*, 227(2), 429–460.
- Demetrescu, M., Rodrigues, P. M. M., & Taylor, A. M. R. (2023).
  Transformed regression-based long-horizon predictability tests.
  *Journal of Econometrics*, 237(2), 105316.
- Harvey, D. I., Leybourne, S. J., & Taylor, A. M. R. (2021). Simple
  tests for stock return predictability with good size and power
  properties. *Journal of Econometrics*, 224(1), 198–214.
- Hosseinkouchack, M., & Demetrescu, M. (2021). Finite-sample size
  control of IVX-based tests in predictive regressions. *Econometric
  Theory*, 37(4), 769–793.
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). Robust
  econometric inference for stock return predictability. *Review of
  Financial Studies*, 28(5), 1506–1553.
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023). Taking
  stock of long-horizon predictability tests: Are factor returns
  predictable? *Journal of Econometrics*, 237(2), 105380.
- Lee, J. H. (2016). Predictive quantile regression with persistent
  covariates: IVX-QR approach. *Journal of Econometrics*, 192(1),
  105–118.
- Magdalinos, T. (2022). Least squares and IVX limit theory in systems
  of predictive regressions with GARCH innovations. *Econometric
  Theory*, 38(5), 875–912.
- Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the
  predictability of US housing price index returns based on an IVX-AR
  model. *Journal of the American Statistical Association*, 115(532),
  1598–1619.

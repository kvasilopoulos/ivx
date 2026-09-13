# IVX predictive regressions

``` r

library(ivx)
```

This vignette covers the baseline method of the package: the IVX
estimator and Wald test of Kostakis, Magdalinos & Stamatogiannis (2015,
*KMS*), its long-horizon version (Kostakis, Magdalinos & Stamatogiannis
2023), and the tuning and standard-error options exposed by
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md).

## The problem

The predictive regression

\\ y_t = \mu + \beta' x\_{t-1} + u_t, \qquad x_t = \mu_x + R_n
x\_{t-1} + v_t, \qquad t = 1, \dots, n, \\

relates a stationary variable (a return) to lagged predictors that are
often highly persistent, \\R_n = I + C / n^{\alpha}\\ with \\\alpha \in
(0, 1\]\\, and whose innovations \\v_t\\ are correlated with \\u_t\\. In
that setting the OLS t-statistic is not asymptotically normal: its
distribution depends on the unknown local-to-unity parameter \\C\\ and
on the endogeneity correlation \\\delta = \mathrm{corr}(u_t, v_t)\\, and
it over-rejects. The package reports \\\delta\\ through
[`delta()`](https://kvasilopoulos.github.io/ivx/reference/delta.md).

## The IVX instrument

KMS instrument \\x\_{t-1}\\ with a self-generated variable of
controlled, mildly integrated persistence,

\\ z_t = \sum\_{j=0}^{t-1} \rho_n^{\\j} \Delta x\_{t-j}, \qquad \rho_n =
1 - \frac{c_z}{n^{\beta}}, \quad \beta \in (0, 1), \\ c_z \> 0, \\

i.e. \\z_t = \rho_n z\_{t-1} + \Delta x_t\\ with \\z_0 = 0\\. Whatever
the persistence of \\x_t\\, \\z_t\\ is mildly integrated, and the IV
estimator

\\ \hat\beta\_{ivx} = \Big( \sum_t z\_{t-1} \tilde x\_{t-1}' \Big)^{-1}
\sum_t z\_{t-1} \tilde y_t, \qquad \tilde x\_{t-1} = x\_{t-1} - \bar x,
\\ \tilde y_t = y_t - \bar y, \\

is asymptotically mixed normal. The Wald statistic for \\H_0: H\beta =
h\\,

\\ W = (H\hat\beta\_{ivx} - h)' \big\[ H \hat Q H' \big\]^{-1}
(H\hat\beta\_{ivx} - h) \\\to\\ \chi^2\_{\mathrm{rank}(H)}, \\

uses the covariance \\\hat Q = (Z'\tilde X)^{-1} M (\tilde X' Z)^{-1}\\
with, in KMS’s finite-sample corrected form,

\\ M = \hat\sigma_u^2 \sum_t z\_{t-1} z\_{t-1}' \\-\\ n \\ \bar z \bar
z' \\ \hat\Sigma\_{FM}, \qquad \hat\Sigma\_{FM} = \hat\sigma_u^2 -
\hat\Omega\_{uv}' \hat\Omega\_{vv}^{-1} \hat\Omega\_{uv}, \\

where \\\hat\Omega\\ are Newey–West long-run (co)variances of the OLS
residuals \\\hat u_t\\ and of the AR(1) residuals of each predictor,
with Bartlett bandwidth \\\lfloor n^{1/3} \rfloor\\. The subtracted term
corrects the finite-sample effect of estimating the intercept (KMS,
p. 1516) and is what the package reports;
[`summary()`](https://rdrr.io/r/base/summary.html) prints the individual
Wald statistics (equal to the squared t-ratios) and the joint statistic.

``` r

mod <- ivx(Ret ~ DP + TBL, data = kms)
summary(mod)
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
delta(mod)
#> [1] -0.9755750 -0.0610432
```

## Tuning

KMS recommend \\\beta = 0.95\\ and \\c_z = 1\\, the package defaults.
Both, and the Newey–West bandwidth, are arguments of
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md):

``` r

coef(ivx(Ret ~ DP + TBL, data = kms, beta = 0.9, cz = 5, bandwidth = 10))
#>          DP         TBL 
#>  0.00501996 -0.12896546
```

Smaller \\\beta\\ (or larger \\c_z\\) makes the instrument less
persistent: better size under strong endogeneity at some cost in power.
Hosseinkouchack & Demetrescu (2021) show the convergence to the
\\\chi^2\\ limit is slower the closer \\\rho_n\\ is to one, and Lee
(2016) uses \\c_z = 5\\ for quantile regressions.

## Long horizons

For a horizon \\K \> 1\\ the regression is between the \\K\\-period sum
\\y_t(K) = \sum\_{i=0}^{K-1} y\_{t+i}\\ and \\x\_{t-1}\\. Overlapping
sums make the usual long-horizon OLS statistics oversized. The package
implements the IVX-Wald statistic of Kostakis, Magdalinos &
Stamatogiannis (2023), eqs (15) and (23): the estimator regresses
\\y_t(K)\\ on the \\K\\-period sum \\x\_{t-1}(K)\\ using the single-lag
instrument \\z\_{t-1}\\,

\\ \tilde A_K = Y(K)' Z \\\big\[ X(K)' Z \big\]^{-1}, \\

while the covariance uses the \\K\\-period sum of the instrument,
\\z\_{t-1}(K)\\, in place of \\z\_{t-1}\\ in \\M\\ above, with
\\\hat\sigma^2_u\\ still from the one-step regression. The statistic
keeps its \\\chi^2\\ limit.

``` r

summary(ivx(Ret ~ DP + TBL, data = kms, horizon = 12))
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 12)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)  
#> DP   0.008202   0.004760   1.723    2.970    0.0848 .
#> TBL -0.062903   0.059869  -1.051    1.104    0.2934  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Joint Wald statistic:  3.998 on 2 DF, p-value 0.1355
#> Multiple R-squared:  0.05429,    Adjusted R-squared:  0.06255
```

## Lag augmentation for power

The IVX instrument is less persistent than \\x_t\\, so the test loses
local power against OLS exactly when \\x_t\\ is near-integrated and
endogenous. Demetrescu (2014) shows that adding \\y\_{t-1}\\ to the
regression, instrumented by itself, \\ y_t = \phi\\ y\_{t-1} + \beta'
x\_{t-1} + u_t, \qquad \phi = 0 \text{ under the null}, \\ feeds the
signal back into the instrument and can raise power substantially when
the instrument is weak (small \\\eta\\), while being asymptotically
equivalent to plain IVX otherwise. `lag_y = TRUE` does this; the joint
Wald statistic still tests only \\\beta\\.

``` r

ivx(Ret ~ DP + TBL, data = kms, lag_y = TRUE)
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, lag_y = TRUE, horizon = 1)
#> 
#> Coefficients:
#>        DP        TBL      y_lag  
#>  0.007248  -0.073288   0.091950
```

## Heteroskedasticity-robust standard errors

Demetrescu, Georgiev, Rodrigues & Taylor (2023) show that the IVX
statistics keep their standard limits under much weaker assumptions on
the innovations (unconditional and conditional heteroskedasticity) if
Eicker–White standard errors are used, i.e. \\\hat\sigma_u^2 \sum z z'\\
is replaced by \\\sum z\_{t-1} z\_{t-1}' \hat u_t^2\\:

``` r

summary(ivx(Ret ~ DP + TBL, data = kms, robust = TRUE))
#> 
#> Call:
#> ivx(formula = Ret ~ DP + TBL, data = kms, robust = TRUE, horizon = 1)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP   0.006145   0.004792   1.282    1.644     0.200
#> TBL -0.080717   0.057250  -1.410    1.988     0.159
#> (Eicker-White standard errors)
#> 
#> Joint Wald statistic:  2.893 on 2 DF, p-value 0.2354
#> Multiple R-squared:  0.004968,   Adjusted R-squared:  0.003036
```

`robust = TRUE` is available for `horizon = 1`; the long-horizon
heteroskedasticity-robust form is not defined in the literature. See
[`vignette("robust-inference")`](https://kvasilopoulos.github.io/ivx/articles/robust-inference.md)
for the bootstrap alternative.

## Replication of Kostakis et al. (2015)

The `kms` dataset is the authors’ monthly file. The published estimates
(Table 6, p. 1531 and Table 8, p. 1537) are reproduced to the printed
precision, and these values are asserted in the package tests.

``` r

m6 <- ivx(Ret ~ DE, data = kms)
c(coef = round(coef(m6), 4), wald = round(m6$Wald_Ind, 3), delta = round(delta(m6), 3))
#> coef.DE wald.DE   delta 
#> -0.0033  0.3930 -0.0670
# paper: -0.0033, 0.393, -0.067

m8 <- ivx(Ret ~ DP + TBL, data = kms)
c(round(coef(m8), 4), joint = round(m8$Wald_Joint, 3))
#>      DP     TBL   joint 
#>  0.0061 -0.0807  3.6440
# paper: 0.0061, -0.0807, 3.644
```

Long-horizon Wald statistics for `Ret ~ EP + TBL` (Table 13, p. 1547):

``` r

h <- c(4, 12, 24, 36, 48, 60)
t13 <- t(sapply(h, function(k) {
  m <- ivx(Ret ~ EP + TBL, data = kms, horizon = k)
  c(EP = m$Wald_Ind[["EP"]], TBL = m$Wald_Ind[["TBL"]], joint = m$Wald_Joint)
}))
round(cbind(horizon = h, t13), 3)
#>      horizon    EP   TBL joint
#> [1,]       4 5.778 3.894 7.638
#> [2,]      12 6.383 3.166 7.614
#> [3,]      24 4.990 2.124 5.794
#> [4,]      36 4.599 1.915 5.383
#> [5,]      48 4.983 1.441 5.660
#> [6,]      60 4.321 1.039 4.822
# paper: EP 5.778 6.383 4.990 4.599 4.983 4.321 | TBL 3.894 3.166 2.124 1.915 1.441 1.039
```

## Caveats

- The regression must carry an intercept;
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)
  handles it internally and warns if the formula removes it.
- The instrument uses the whole sample, so
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) on a
  subsample is not the same as a subsample statistic with the
  full-sample instrument; use
  [`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md)
  for the latter.
- With strong endogeneity (\\\|\delta\|\\ near 1) and near-unit-root
  predictors the asymptotic test still over-rejects in samples of a few
  hundred observations, especially one-sided.
  [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  and
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  address this.
- Weighted fits (`weights`) are supported but rarely appropriate for the
  theory.

## References

- Demetrescu, M. (2014). Enhancing the local power of IVX-based tests in
  predictive regressions. *Economics Letters*, 124(2), 269–273.

- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). Robust
  econometric inference for stock return predictability. *Review of
  Financial Studies*, 28(5), 1506–1553.

- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023). Taking
  stock of long-horizon predictability tests: Are factor returns
  predictable? *Journal of Econometrics*, 237(2), 105380.

- Magdalinos, T., & Phillips, P. C. B. (2009). Limit theory for
  cointegrated systems with moderately integrated and moderately
  explosive regressors. *Econometric Theory*, 25(2), 482–526.

- Phillips, P. C. B., & Lee, J. H. (2013). Predictive regression under
  various degrees of persistence and robust long-horizon regression.
  *Journal of Econometrics*, 177(2), 250–264.

- Hosseinkouchack, M., & Demetrescu, M. (2021). Finite-sample size
  control of IVX-based tests in predictive regressions. *Econometric
  Theory*, 37(4), 769–793.

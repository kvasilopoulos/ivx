# Non-IVX benchmarks: Bonferroni Q, ARM, hybrid t, empirical likelihood

``` r

library(ivx)
```

Five non-IVX procedures are included as the benchmarks the IVX
literature compares against: the Bonferroni Q-test of Campbell & Yogo
(2006),
[`cy_test()`](https://kvasilopoulos.github.io/ivx/reference/cy_test.md);
the augmented regression method of Amihud, Hurvich & Wang (2009),
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md); the
control-function regression of Elliott (2011),
[`elliott_cf()`](https://kvasilopoulos.github.io/ivx/reference/elliott_cf.md);
the hybrid switching t-test of Harvey, Leybourne & Taylor (2021),
[`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md);
and the unified empirical likelihood test of Liu, Yang, Cai & Peng
(2019),
[`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md).

## Bonferroni Q-test

[`cy_test()`](https://kvasilopoulos.github.io/ivx/reference/cy_test.md)
is the Campbell & Yogo (2006) procedure, the feasible version of the
sup-bound / Bonferroni idea of Cavanagh, Elliott & Stock (1995). If
\\\rho\\ were known, the UMP conditional test of \\\beta = 0\\ is the
t-ratio of the regression augmented with \\x_t - \rho x\_{t-1}\\ (the
“Q-test”), whose estimate \\ \hat\beta(\rho) = \frac{\sum_t
x^\mu\_{t-1}\big(y_t - \tfrac{\sigma\_{ue}}{\sigma_e\omega}(x_t - \rho
x\_{t-1})\big) -
\tfrac{T}{2}\tfrac{\sigma\_{ue}}{\sigma_e\omega}(\omega^2 -
\sigma_v^2)}{\sum_t x^{\mu 2}\_{t-1}} \\ (their eq. 25, with the
AR(\\p\\) correction of Appendix A) has standard error
\\\sigma_u\sqrt{1-\delta^2}/(\sum x^{\mu 2}\_{t-1})^{1/2}\\. Since
\\\rho\\ is not consistently estimable, a confidence interval for it is
obtained by inverting the DF-GLS statistic (Stock, 1991; Elliott,
Rothenberg & Stock, 1996) and the Bonferroni interval for \\\beta\\ runs
from \\\hat\beta(\bar\rho) - 1.645\\se\\ to
\\\hat\beta(\underline\rho) + 1.645\\se\\. The levels of the interval
for \\\rho\\ come from the paper’s Table 2, which tightens the plain
Bonferroni bound so that the one-sided test has size exactly 5% for some
\\c\\.

The DF-GLS null quantiles as a function of \\c\\ are simulated once
(`data-raw/dfgls-quantiles.R`, \\c \in \[-100, 10\]\\, 20 000
replications of a 600-step OU process) and shipped as internal data; the
5% quantile at \\c = 0\\ reproduces the \\-1.95\\ of Elliott et
al. (1996).

``` r

cy_test(Ret ~ DP, data = kms)
#> 
#> Call:
#> cy_test(formula = Ret ~ DP, data = kms)
#> 
#> Bonferroni Q-test (Campbell & Yogo, 2006)
#> 
#> delta = -0.972, DF-GLS = -1.468 (p = 2), CI for c at levels (0.055, 0.082): [-9.319, 1.044], rho: [0.991, 1.001]
#> OLS slope = 0.006128; Q-estimates at the ends of the rho interval: 0.009069, 0.0004869
#> 90% Bonferroni confidence interval for beta: [-0.0009783, 0.01053]
#> 5% one-sided Q-tests: H1 beta > 0 do not reject H0; H1 beta < 0 do not reject H0
cy_test(Ret ~ EP, data = kms)
#> 
#> Call:
#> cy_test(formula = Ret ~ EP, data = kms)
#> 
#> Bonferroni Q-test (Campbell & Yogo, 2006)
#> 
#> delta = -0.7909, DF-GLS = -3.014 (p = 3), CI for c at levels (0.065, 0.17): [-29.24, -11.74], rho: [0.9717, 0.9886]
#> OLS slope = 0.008698; Q-estimates at the ends of the rho interval: 0.02199, 0.0147
#> 90% Bonferroni confidence interval for beta: [0.01057, 0.02612]
#> 5% one-sided Q-tests: H1 beta > 0 reject H0; H1 beta < 0 do not reject H0
```

In a Monte Carlo with \\T = 200\\ and 500 replications, the right-tailed
5% Q-test rejects a true null 4.4–5.8% of the time and the left-tailed
one 1.2–3.6% for \\\delta \in \\-0.9, -0.5\\\\ and \\c \in \\0, -20\\\\
— the asymmetry the paper describes (Section 3.4: left-tailed
probability “can be as small as 1.2%”).

## Augmented regression method

[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) is
reduced-bias OLS for *stationary* persistent predictors (Demetrescu &
Rodrigues, 2022, build
[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md) on
the same augmentation idea).

### Idea

With \\y_t = \alpha + \beta' x\_{t-1} + u_t\\ and \\x_t = \mu + \Phi
x\_{t-1} + v_t\\, write \\u_t = \phi' v_t + e_t\\ with \\e_t\\
independent of \\x\\. OLS of \\y_t\\ on \\x\_{t-1}\\ is biased because
\\u_t\\ is correlated with \\v_t\\ and hence with the future of the
regressor (Stambaugh, 1999). Regressing \\y_t\\ on \\x\_{t-1}\\ *and*
\\v_t\\ removes it — if \\v_t\\ were known. Using OLS VAR residuals does
not help (they are orthogonal to \\x\_{t-1}\\ by construction), so the
VAR coefficient matrix is first bias-corrected with the Nicholls & Pope
(1988) expansion \\ E\[\hat\Phi - \Phi\] =
-\tfrac{1}{n}\\\Sigma_v\Big\[(I-\Phi')^{-1} + \Phi'(I-\Phi'^2)^{-1} +
\sum\_{\lambda \in \mathrm{Spec}(\Phi)} \lambda
(I-\lambda\Phi')^{-1}\Big\]\Gamma_0^{-1}, \\ iterated up to `iter` times
(Kendall’s \\1 + 3\phi\\ in the scalar case), and the corrected
residuals \\v^c_t = x_t - \hat\mu^c - \hat\Phi^c x\_{t-1}\\ are the
control variables.

The OLS covariance of the augmented regression is wrong for
\\\hat\beta^c\\ because \\v^c_t\\ is estimated. The paper’s estimator
(eqs 7–8) adds the uncertainty of \\\hat\Phi^c\\: \\
\widehat{\mathrm{cov}}\[\hat\beta^c\] =
(\hat\phi'\hat\Sigma_v\hat\phi)\\(X'X)^{-1}\_{xx} + \hat\sigma_e^2
\Big\[\tfrac{\sum_t r\_{it} r\_{jt}}{\sum_t r\_{it}^2 \sum_t
r\_{jt}^2}\Big\]\_{ij}, \\ with \\r\_{jt}\\ the residual of
\\x\_{j,t-1}\\ on the other regressors of the augmented regression.
[`summary()`](https://rdrr.io/r/base/summary.html) reports \\t\\-ratios
and the joint Wald test from this matrix.

``` r

m <- arm(Ret ~ DP + TBL, data = kms)
summary(m)
#> 
#> Call:
#> arm(formula = Ret ~ DP + TBL, data = kms)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP   0.002050   0.003797   0.540    0.292     0.589
#> TBL -0.046509   0.055999  -0.831    0.690     0.406
#> 
#> Joint Wald statistic:  1.068 on 2 DF, p-value 0.5862
#> Multiple R-squared:  0.02415,    Adjusted R-squared:  0.0194
m$Phi      # bias-corrected VAR(1) coefficients
#>                DP           TBL
#> DP   9.965559e-01 -0.0001402752
#> TBL -5.682961e-05  0.9964614118
m$Phi_ols  # OLS
#>               DP        TBL
#> DP   0.992685362 0.02492259
#> TBL -0.000183203 0.99299626
```

### Check

With \\\phi = 0.9\\, \\\mathrm{corr}(u, v) = -0.9\\, \\n = 100\\ and
\\\beta = 0\\ the OLS slope has mean bias about \\0.037\\ (\\\approx
-\rho(1+3\phi)/n\\) over 400 replications;
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) reduces
it to \\0.003\\, with \\\hat\Phi^c\\ averaging \\0.894\\ against an OLS
mean of \\0.857\\. In the paper’s two-predictor Case 1 the 5%
\\t\\-tests reject 7–8% at \\n = 100\\, in line with the 6–10% (\\n =
50\\) and 5–7% (\\n = 200\\) the paper reports.

### Caveats

- Assumes stationary predictors: the correction diverges as the largest
  eigenvalue of \\\Phi\\ approaches one, and the iteration stops if the
  corrected VAR becomes non-stationary. For near-unit-root predictors
  use [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) or
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md).
- VAR(1) only, short horizon only, no bootstrap.

## Control-function regression

Elliott (2011) shows that if stationary covariates \\z_t\\ are available
that are contemporaneously correlated with the shocks to both the
predictor and the response, the Wald test of \\eta = 0\\ in \\ y_t =
lpha + eta' x\_{t-1} + \gamma' Z_t + ilde u_t, \qquad Z_t = (z_t',
z\_{t-1}', \dots, z\_{t-q}')', \\ is asymptotically \\\chi^2\\ whatever
the persistence of \\x_t\\, because the covariates “orthogonalise” the
innovations (his Theorem 2 versus the Elliott–Stock 1994 distribution of
Theorem 1). The covariates are a modelling choice, not a data
construction —
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) and
[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
are the feasible versions where the control variable is built from the
predictor’s own innovations.
[`elliott_cf()`](https://kvasilopoulos.github.io/ivx/reference/elliott_cf.md)
runs the regression, reports the Wald test with Eicker-White standard
errors and the correlation that remains between the regression residuals
and the predictor innovations, which should be near zero.

``` r

elliott_cf(Ret ~ DP, ~ TBL, data = kms, lags = 1)
#> 
#> Call:
#> elliott_cf(formula = Ret ~ DP, covariates = ~TBL, data = kms, 
#>     lags = 1)
#> 
#> Control-function predictive regression (Elliott, 2011), 1 covariate lag(s)
#> 
#>    Estimate Std. Error t value Pr(>|t|)
#> DP 0.005605   0.005102   1.099    0.272
#> (Eicker-White standard errors)
#> 
#> Wald statistic: 1.207 on 1 DF, p-value 0.2719
#> Remaining innovation correlation: DP -0.977
```

With a covariate that does absorb the correlation (\$z_t = v_t + \$
noise, \\\delta = -0.9\\, \\c = 0\\, \\T = 200\\) the 5% test rejects a
true null 4% of the time against 23% for plain OLS (500 replications).

## Hybrid t-test

[`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md)
implements \\T\_{hyb}\\ of Harvey, Leybourne & Taylor (2021) for a
single predictor. It uses only regression t-ratios: the standard one,
\\T\\, and \$ ilde T\\, in which the predictor is quasi-GLS demeaned
with\\ar c = 7\$ as in Elliott et al. (1996). Under strong persistence
their null distributions depend on \\c\\ and on \$ ho\_{xy}\$, so the
paper tabulates *conservative* critical values — the maximum over \\c\\
of the \\lpha\\-level quantile — as response surfaces in \\\hat
ho\_{xy}\\ (their Table 1), and switches:

1.  if the ADF normalised-bias statistic (lag length by MBIC) is below
    \\-4\sqrt{T}\\ the predictor is weakly persistent: \\T\\ with a
    normal critical value;
2.  otherwise, for an upper-tail test, \\T\\ with \\cv(\hat ho\_{xy})\\
    when \\\hat ho\_{xy} \> -0.1\\ and \$ ilde T\$ with \$ ilde{cv}(
    ho\_{xy})\$ when \\\hat ho\_{xy} \< -0.1\\ (quasi-GLS demeaning pays
    off only when the endogeneity is strong); lower-tail tests mirror
    this.

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
hlt_test(Ret ~ TBL, data = kms, alternative = "less", level = 0.1)
#> 
#> Call:
#> hlt_test(formula = Ret ~ TBL, data = kms, alternative = "less", 
#>     level = 0.1)
#> 
#> Hybrid predictability test (Harvey, Leybourne & Taylor, 2021)
#> 
#> Selected test: T_con (OLS-demeaned t, conservative critical value)
#> statistic = -1.403, 10% critical value = -1.306 (alternative: beta < 0): reject
#> slope = -0.07836, t = -1.403, t (quasi-GLS) = -1.384, ADF = -6.707 (p = 6, cutoff -128.6), rho_xy = -0.05541
```

In the paper’s design (\\T = 200\\, upper-tail 5% test) the empirical
size in 1000 replications is 0.07/0.04/0.05 for \$ ho\_{xy} = -0.9\$ and
\$ ho = 1, 0.95, 0.5\$, and 0.00/0.02/0.04 for \$ ho\_{xy} = 0.9\$ — the
“undersize for positive \$ ho\_{xy}\$ in the more persistent cases” the
paper notes is the price of the conservative critical values. The test
is one-sided by construction and has no p-value; the object reports the
selected test, its statistic and critical value.

## Unified empirical likelihood test

[`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md)
implements Liu, Yang, Cai & Peng (2019). Their model adds the lagged
difference of the predictor, \\ Y_t = lpha + eta_1 \Delta X\_{t-1} +
eta_2 X\_{t-2} + U_t, \\ so that \\Y_t\\ can be stationary whether or
not \\X_t\\ is. The intercept is removed by differencing at lag \\m =
\lfloor n/2 floor\\ (Zhu, Cai & Peng, 2014), and the empirical
likelihood is built on the two score equations, the second weighted by
\\1/\sqrt{1 + ilde X\_{t-2}^2}\\ so that its sample variance converges
whatever the persistence. The profile EL ratios for \\eta_2 = 0\\ (no
predictability), \\eta_1 = 0\\ and the joint null are then
\\\chi^2(1)\\, \\\chi^2(1)\\, \\\chi^2(2)\\ with no tuning parameter and
no persistence classification (Theorem 2). The price is efficiency: only
\\m - 2\\ differenced observations enter, and the paper reports that the
test on \\eta_1\\ is oversized at \\n = 200\\.

``` r

el_test(Ret ~ DP, data = kms)
#> 
#> Call:
#> el_test(formula = Ret ~ DP, data = kms)
#> 
#> Unified empirical likelihood test (Liu, Yang, Cai & Peng, 2019), m = 516
#> 
#>                  Estimate EL ratio    df Pr(> chi)
#> beta1 (dX[t-1]) -0.133221    4.456 1.000    0.0348
#> beta2 (X[t-2])   0.008196    2.404 1.000    0.1210
#> joint                        7.250 2.000    0.0266
```

The EL dual problem is solved by Newton’s method with Owen’s (2001)
pseudo-logarithm; the profiles are minimised over the nuisance
coefficient by BFGS from the OLS start.

### References

- Amihud, Y., Hurvich, C. M., & Wang, Y. (2009). Multiple-predictor
  regressions: Hypothesis testing. *Review of Financial Studies*, 22(1),
  413–434.
- Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
  reduced-bias estimation method. *Journal of Financial and Quantitative
  Analysis*, 39(4), 813–841.
- Campbell, J. Y., & Yogo, M. (2006). Efficient tests of stock return
  predictability. *Journal of Financial Economics*, 81(1), 27–60.
- Cavanagh, C. L., Elliott, G., & Stock, J. H. (1995). Inference in
  models with nearly integrated regressors. *Econometric Theory*, 11(5),
  1131–1147.
- Elliott, G. (2011). A control function approach for testing the
  usefulness of trending variables in predictive regressions and
  econometric models. *Journal of Econometrics*, 164(1), 79–91.
- Elliott, G., Rothenberg, T. J., & Stock, J. H. (1996). Efficient tests
  for an autoregressive unit root. *Econometrica*, 64(4), 813–836.
- Harvey, D. I., Leybourne, S. J., & Taylor, A. M. R. (2021). Simple
  tests for stock return predictability with good size and power
  properties. *Journal of Econometrics*, 224(1), 198–214.
- Liu, X., Yang, B., Cai, Z., & Peng, L. (2019). A unified test for
  predictability of asset returns regardless of properties of predicting
  variables. *Journal of Econometrics*, 208(1), 141–159.
- Nicholls, D. F., & Pope, A. L. (1988). Bias in the estimation of
  multivariate autoregressions. *Australian Journal of Statistics*, 30A,
  296–309.
- Owen, A. B. (2001). *Empirical Likelihood*. Chapman & Hall.
- Stock, J. H. (1991). Confidence intervals for the largest
  autoregressive root in U.S. macroeconomic time series. *Journal of
  Monetary Economics*, 28(3), 435–459.
- Stambaugh, R. F. (1999). Predictive regressions. *Journal of Financial
  Economics*, 54(3), 375–421.

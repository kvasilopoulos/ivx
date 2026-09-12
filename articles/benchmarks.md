# Non-IVX benchmarks: ARM, hybrid t-test, empirical likelihood

``` r

library(ivx)
```

Three non-IVX procedures are included as the benchmarks the IVX
literature compares against: the augmented regression method of Amihud,
Hurvich & Wang (2009),
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md), the
hybrid switching t-test of Harvey, Leybourne & Taylor (2021),
[`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md),
and the unified empirical likelihood test of Liu, Yang, Cai & Peng
(2019),
[`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md).

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
- Stambaugh, R. F. (1999). Predictive regressions. *Journal of Financial
  Economics*, 54(3), 375–421.

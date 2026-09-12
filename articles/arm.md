# Augmented regression method (non-IVX benchmark)

``` r

library(ivx)
```

[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md)
implements the multipredictor augmented regression method (mARM) of
Amihud, Hurvich & Wang (2009). It is not an IVX estimator: it is
reduced-bias OLS for *stationary* persistent predictors, and is included
as the benchmark the IVX literature compares against (Demetrescu &
Rodrigues, 2022, build
[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md) on
the same augmentation idea).

## Idea

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

## Check

With \\\phi = 0.9\\, \\\mathrm{corr}(u, v) = -0.9\\, \\n = 100\\ and
\\\beta = 0\\ the OLS slope has mean bias about \\0.037\\ (\\\approx
-\rho(1+3\phi)/n\\) over 400 replications;
[`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) reduces
it to \\0.003\\, with \\\hat\Phi^c\\ averaging \\0.894\\ against an OLS
mean of \\0.857\\. In the paper’s two-predictor Case 1 the 5%
\\t\\-tests reject 7–8% at \\n = 100\\, in line with the 6–10% (\\n =
50\\) and 5–7% (\\n = 200\\) the paper reports.

## Caveats

- Assumes stationary predictors: the correction diverges as the largest
  eigenvalue of \\\Phi\\ approaches one, and the iteration stops if the
  corrected VAR becomes non-stationary. For near-unit-root predictors
  use [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) or
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md).
- VAR(1) only, short horizon only, no bootstrap.

## References

- Amihud, Y., Hurvich, C. M., & Wang, Y. (2009). Multiple-predictor
  regressions: Hypothesis testing. *Review of Financial Studies*, 22(1),
  413–434.
- Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
  reduced-bias estimation method. *Journal of Financial and Quantitative
  Analysis*, 39(4), 813–841.
- Nicholls, D. F., & Pope, A. L. (1988). Bias in the estimation of
  multivariate autoregressions. *Australian Journal of Statistics*, 30A,
  296–309.
- Stambaugh, R. F. (1999). Predictive regressions. *Journal of Financial
  Economics*, 54(3), 375–421.

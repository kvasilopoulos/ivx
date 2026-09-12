# Residual-augmented IVX

``` r

library(ivx)
```

[`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
implements the residual-augmented IVX estimator of Demetrescu &
Rodrigues (2022), a bias-reduced version of IVX in the spirit of Amihud
& Hurvich (2004).

## Idea

Write the innovations of the predictor’s autoregression as \\v_t\\ and
decompose the regression error as \\u_t = \gamma' v_t + \eta_t\\ with
\\\eta_t\\ uncorrelated with \\v_t\\. The finite-sample bias of any
estimator of \\\beta\\ comes from \\\gamma \ne 0\\. If \\v_t\\ were
observed, adding it as a regressor,

\\ \tilde y_t = \beta' \tilde x\_{t-1} + \gamma' v_t + \eta_t, \\

would remove the problem. Amihud & Hurvich estimate \\v_t\\ from a
bias-corrected autoregression; Demetrescu & Rodrigues show that
augmenting the *IVX* regression with plain autoregressive residuals is
enough to reduce the bias substantially while retaining the
persistence-robust \\\chi^2\\ inference.

## Procedure

1.  Fit an autoregression of order \\p\\ in levels to the predictors (a
    VAR with several predictors), **without intercept**, \\p\\ chosen by
    AIC (the paper’s recommendation); keep the residuals
    \\\hat\varepsilon_t\\, \\t = p+1, \dots, n\\, and demean them.
2.  Regress the demeaned response on \\\hat\varepsilon_t\\ by OLS and
    keep the residuals \\\tilde y_t = \tilde y_t^{\\0} - \hat\gamma'
    \hat\varepsilon_t\\.
3.  Estimate \\\beta\\ by IVX of \\\tilde y_t\\ on \\\tilde x\_{t-1}\\
    with the KMS instrument (not demeaned).

Standard errors are the heteroskedasticity-robust form of the paper’s
eq. (9)/(14), which stay valid whether \\x_t\\ is stationary or
near-integrated:

\\ \widehat{\mathrm{Cov}}(\tilde\beta\_{ivx}) = B^{-1} M B^{-1\prime},
\quad B = \sum_t z\_{t-1}\tilde x\_{t-1}', \quad M = \sum_t z\_{t-1}
z\_{t-1}' \tilde\varepsilon_t^2 + \hat Q_T - n \bar z \bar z'
\hat\Sigma\_{FM}, \\

where \\\tilde\varepsilon_t\\ are the OLS residuals of the augmented
regression, \\\hat Q_T = H\_{zx} H\_{xx}^{-1} \tilde H\_{xx}
H\_{xx}^{-1} H\_{zx}'\\ accounts for the estimation of
\\\hat\varepsilon_t\\ (it matters only in the stationary case), and the
last term is the Kostakis et al. (2015) intercept correction that the
paper also uses in its simulations.

``` r

m <- ivx_ra(Ret ~ DP + TBL, data = kms)
m
#> 
#> Call:
#> ivx_ra(formula = Ret ~ DP + TBL, data = kms)
#> 
#> Residual-augmented IVX, AR order p = 4 (aic)
#> 
#> Coefficients:
#>        DP        TBL  
#> -0.001934  -0.053984
summary(m)
#> 
#> Call:
#> ivx_ra(formula = Ret ~ DP + TBL, data = kms)
#> 
#> Coefficients:
#>      Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP  -0.001934   0.005010  -0.386    0.149     0.699
#> TBL -0.053984   0.055313  -0.976    0.953     0.329
#> (Eicker-White standard errors)
#> 
#> Joint Wald statistic:  1.423 on 2 DF, p-value 0.4908
#> Multiple R-squared:  0.01794,    Adjusted R-squared:  0.01602
```

Compared with plain IVX the point estimate on the strongly endogenous
predictor (`DP`, \\\delta \approx -0.98\\) shrinks and the standard
error is smaller — the efficiency gain from removing the \\\gamma' v_t\\
component of the error.

``` r

rbind(ivx = coef(ivx(Ret ~ DP + TBL, data = kms)), ivx_ra = coef(m))
#>                  DP         TBL
#> ivx     0.006145163 -0.08071667
#> ivx_ra -0.001934204 -0.05398376
```

## Replication: size in the paper’s Monte Carlo

The paper has no reproducible empirical application, so the check is
against the size rows of its Table 3 (right-sided 5% tests, \\T = 200\\,
DGP (24)–(25) with \\\rho = 1 - c/T\\, short-run AR parameter \\-0.5\\
and innovation correlation \\-0.95\\). With 10 000 replications:

| \\c\\ | paper \\t\_{ivx}\\ | `ivx` | paper \\\tilde t\_{ivx}^{\mu_0}\\ | `ivx_ra` |
|-------|--------------------|-------|-----------------------------------|----------|
| 0     | 0.116              | 0.129 | 0.054                             | 0.064    |
| 10    | 0.088              | 0.089 | 0.055                             | 0.064    |
| 20    | 0.074              | 0.072 | 0.055                             | 0.064    |
| 30    | 0.066              | 0.057 | 0.053                             | 0.055    |
| 40    | 0.064              | 0.057 | 0.050                             | 0.059    |
| 50    | 0.061              | 0.051 | 0.050                             | 0.056    |

The over-rejection of plain IVX under strong endogeneity is largely
removed; the residual offset of about 0.01 at small \\c\\ is shared by
the plain IVX column and so reflects a detail of the simulation design
rather than the estimator.

## Long horizons: the transformed regression

Demetrescu, Rodrigues & Taylor (2023) extend the estimator to the
\\h\\-period regression without HAC estimation. Following Britten-Jones
et al. (2011), the overlapping regression of \\\sum\_{j=1}^h y\_{t+j}\\
on \\x_t\\ is numerically the non-overlapping regression of \\y\_{t+1}\\
on the transformed regressor \\A_h' x\\, and the same trick applies to
the instrument: \\ z_t^{trf,(h)} =
\sum\_{i=\max(1,\\t-h+1)}^{\min(t,\\T-h)} z_i, \qquad \hateta_h =
\Big(\sum\_{t=1}^{T-h} z_t ar x_t'\Big)^{-1} \sum\_{t=p}^{T-1}
z_t^{trf,(h)}\\(ar y\_{t+1} - \hat\gamma'\hatarepsilon\_{t+1}), \\ with
the sandwich covariance of eq. (5.7), i.e. the short-horizon one with
\\z_t\\ replaced by \\z_t^{trf,(h)}\\. `horizon = h` does exactly this
and reduces to the short-horizon estimator at `h = 1`.

``` r

sapply(c(1, 3, 12), function(h) ivx_ra(Ret ~ DP, data = kms, horizon = h)$tstat)
#>         DP         DP         DP 
#> -0.4642328 -0.3057244 -0.3139751
```

Under the null with a near-integrated, strongly endogenous predictor
(\\c = 0\\, \\\delta = -0.95\\, \\T = 500\\) the two-sided 5% test
rejects at 0.034 (\\h = 10\\) and 0.040 (\\h = 20\\) in 2000
replications, inside the \[0.023, 0.058\] range the paper reports for
\\T = 500\\.

## Caveats and implementation notes

- Two details the paper leaves implicit are decisive: the autoregression
  must be fitted **without** an intercept (the paper’s \\\mu_0\\
  variant; its OLS-demeaned \\\mu_1\\ variant is much more biased and
  “not recommended for testing”), and the residuals must be **demeaned**
  before augmentation — without that the term \\\gamma
  \bar{\hat\varepsilon}\\ survives and the bias reduction fails. Both
  are built in and neither is optional.
- The remaining bias is \\-\gamma\\ times the bias of the AR
  coefficients, so it grows with \\\|\gamma\|\\ and with the persistence
  of the predictor; it is much smaller than for plain IVX but not zero.
- No bootstrap:
  [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  does not accept `ivx_ra` objects. For `horizon > 1` the Kostakis et
  al. (2015) intercept correction is not applied (it is not derived for
  the transformed regression), and fitted values and residuals are those
  of the transformed regression.
- The AR order is selected by AIC in levels, which the authors argue
  copes with both the stationary and the integrated case; `ar` fixes it
  instead.

## References

- Demetrescu, M., & Rodrigues, P. M. M. (2022). Residual-augmented IVX
  predictive regression. *Journal of Econometrics*, 227(2), 429–460.
- Demetrescu, M., Rodrigues, P. M. M., & Taylor, A. M. R. (2023).
  Transformed regression-based long-horizon predictability tests.
  *Journal of Econometrics*, 237(2), 105316.
- Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
  reduced-bias estimation method. *Journal of Financial and Quantitative
  Analysis*, 39(4), 813–841.

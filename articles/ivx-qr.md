# IVX quantile predictive regression

``` r

library(ivx)
```

[`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
implements the IVX-QR predictability test of Lee (2016): a quantile
regression version of the IVX test that is valid for stationary, mildly
integrated, near-unit-root and mildly explosive predictors alike. It
needs the **quantreg** package.

## Setting

The conditional \\\tau\\-quantile of the return is modelled as

\\ Q\_{y_t}(\tau \mid x\_{t-1}) = \beta\_{0,\tau} + \beta\_\tau'
x\_{t-1}. \\

Ordinary quantile regression on \\x\_{t-1}\\ suffers from the same
problem as OLS: with persistent, endogenous predictors the QR t-ratio
has a non-standard limit that depends on the local-to-unity parameter
and on the quantile-specific endogeneity \\\rho(\tau) =
-\mathrm{corr}\big(1\\u\_{0t} \< 0\\, u\_{xt}\big)\\.

## The test

Lee’s practical procedure (Section 3.3, Proposition 3.2) exploits that
the IVX instrument \\z\_{t-1}\\ is “close” to \\x\_{t-1}\\: run the
quantile regression of \\y_t\\ on an intercept and the demeaned
instrument,

\\ \hat\beta^{IVXQR}\_\tau = \arg\min\_{\beta_0, \beta} \sum_t
\rho\_\tau\big( y_t - \beta_0 - \beta' \tilde z\_{t-1} \big), \\

and use the self-normalised statistic

\\ \frac{\hat f_u(0)^2}{\tau(1-\tau)} \\ \hat\beta\_\tau' \Big(\sum_t
\tilde z\_{t-1} \tilde z\_{t-1}'\Big) \hat\beta\_\tau \\\to\\ \chi^2_K
\quad \text{under } H_0: \beta\_\tau = 0, \\

where \\\hat f_u(0)\\ is the density of the QR residuals at zero,
estimated with a Gaussian kernel and Silverman’s bandwidth (the paper’s
footnote 4). The implied covariance \\\tau(1-\tau)\hat
f_u(0)^{-2}(\tilde Z'\tilde Z)^{-1}\\ gives the standard errors in
[`summary()`](https://rdrr.io/r/base/summary.html).

``` r

m <- ivx_qr(Ret ~ DP, data = kms, tau = 0.5)
summary(m)
#> 
#> Call:
#> ivx_qr(formula = Ret ~ DP, data = kms, tau = 0.5)
#> 
#> IVX-QR at tau = 0.5
#> 
#> Coefficients:
#>    Estimate Std. Error t value Wald Ind Pr(> chi)
#> DP 0.003571   0.005516   0.647    0.419     0.517
#> 
#> Joint Wald statistic:  0.4191 on 1 DF, p-value 0.5174
#> QR endogeneity rho(tau): DP -0.669
```

Several quantiles at once return a list:

``` r

ms <- ivx_qr(Ret ~ DP + TBL, data = kms, tau = c(0.1, 0.25, 0.5, 0.75, 0.9))
t(sapply(ms, function(f) c(tau = f$tau, coef(f), Wald = round(f$Wald_Joint, 2),
                            p = round(1 - pchisq(f$Wald_Joint, 2), 3))))
#>       tau           DP         TBL  Wald     p
#> 0.1  0.10 -0.028305755  0.18938487  6.67 0.036
#> 0.25 0.25 -0.001707079 -0.06790951  0.75 0.688
#> 0.5  0.50  0.008204012 -0.17132111  7.54 0.023
#> 0.75 0.75  0.024622728 -0.19952197 21.90 0.000
#> 0.9  0.90  0.026739798 -0.43020894 34.59 0.000
```

Predictability that is absent at the median can appear in the tails,
which is the empirical point of the paper.

## Tuning

Lee normalises \\c_z = 5\\ (the package default for
[`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md);
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) uses 1)
and chooses the exponent \\\beta\\ from a look-up table indexed by the
estimated endogeneity \\\hat\rho(\tau)\\, so that the asymptotic size of
the nominal 5% test stays below 7.5%: the larger \\\|\hat\rho(\tau)\|\\,
the smaller \\\beta\\. The table is in the paper’s supplement rather
than the paper, so the package reports \\\hat\rho(\tau)\\ (`rho_tau`,
also printed by [`summary()`](https://rdrr.io/r/base/summary.html)) and
leaves `beta` to the user, defaulting to the Kostakis et al. (2015)
value 0.95. In a Monte Carlo with a unit-root predictor and \\\rho =
-0.95\\ the median test with the default has empirical size of about 8%
at \\n = 250\\, which is exactly the situation the rule is meant to
correct; with \\\beta = 0.8\\:

``` r

ivx_qr(Ret ~ DP, data = kms, tau = 0.5, beta = 0.8)$Wald_Joint
#> [1] 0.005426976
```

## Caveats

- Tail quantiles need long samples for an accurate density estimate; the
  paper uses \\n = 700\\ to study the 5% quantile.
- The test in [`summary()`](https://rdrr.io/r/base/summary.html) is the
  QR-on-instrument test of Proposition 3.2, not the full IVX-QR
  estimator of the paper’s equation (3.6), whose non-convex objective
  the paper itself avoids for testing \\\beta\_\tau = 0\\. The
  coefficients are therefore those of the regression on \\\tilde
  z\_{t-1}\\.
- Short horizon only; no bootstrap.
- The `rq` fit is stored in the result (`$rq`) for further quantreg
  methods.

## References

- Lee, J. H. (2016). Predictive quantile regression with persistent
  covariates: IVX-QR approach. *Journal of Econometrics*, 192(1),
  105–118.
- Koenker, R. (2005). *Quantile Regression*. Cambridge University Press.

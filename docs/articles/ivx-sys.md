# Systems of predictive regressions

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ivx`](https://kvasilopoulos.github.io/ivx/)`)`

[`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
estimates several predictive regressions jointly — for instance the
returns on a set of pricing factors on the same predictors, the
empirical setting of Kostakis, Magdalinos & Stamatogiannis (2023) — and
tests restrictions across equations with a single IVX-Wald statistic.

## Model and statistic

\\ y_t = \mu + A\\ x\_{t-1} + \varepsilon_t, \qquad y_t \in
\mathbb{R}^m, \\ x_t \in \mathbb{R}^r, \\

with the predictors as in
[`vignette("ivx")`](https://kvasilopoulos.github.io/ivx/articles/ivx.md)
and \\\varepsilon_t\\ a vector martingale difference with covariance
\\\Sigma\\. The IVX estimator of the \\m \times r\\ matrix \\A\\ (KMS
2023, eq. 15) is

\\ \tilde A_K = Y(K)' Z \\\big\[ X(K)' Z \big\]^{-1}, \\

and the covariance of \\\mathrm{vec}(\tilde A_K)\\ has the Kronecker
form of eq. (23),

\\ \tilde Q_K = \big\[ (Z'X(K))^{-1} \otimes I_m \big\]\\ M_K\\ \big\[
(X(K)'Z)^{-1} \otimes I_m \big\], \qquad M_K = Z(K)'Z(K) \otimes
\hat\Sigma \\-\\ n_K\\ \bar z(K)\bar z(K)' \otimes \hat\Sigma\_{FM}, \\

where \\\hat\Sigma\\ is the covariance of the one-step OLS residuals and
\\\hat\Sigma\_{FM} = \hat\Sigma - \hat\Omega\_{\varepsilon v}'
\hat\Omega\_{vv}^{-1} \hat\Omega\_{\varepsilon v}\\ its
long-run-corrected analogue. For \\H \mathrm{vec}(A) = h\\,

\\ W = \big(H\mathrm{vec}(\tilde A_K) - h\big)' \big\[H \tilde Q_K
H'\big\]^{-1} \big(H\mathrm{vec}(\tilde A_K) - h\big) \\\to\\
\chi^2\_{\mathrm{rank}(H)} . \\

[`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
reports the joint test of \\A = 0\\ (\\mr\\ degrees of freedom), one
test per equation (\\r\\ degrees of freedom) and one per coefficient.
For \\m = 1\\ every number equals the output of
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md), at any
horizon.

\
`s`` ``<-`` `[`ivx_sys`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(``Ret``, ``DE``)`` ``~`` ``DP`` ``+`` ``TBL``, data ``=`` ``kms``)`\
`s`\
`#> `\
`#> Call:`\
`#> ivx_sys(formula = cbind(Ret, DE) ~ DP + TBL, data = kms, horizon = 1)`\
`#> `\
`#> Coefficients (responses in rows):`\
`#>      DP         TBL      `\
`#> Ret   0.006145  -0.080717`\
`#> DE    0.311083  -4.014921`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``s``)`\
`#> `\
`#> Call:`\
`#> ivx_sys(formula = cbind(Ret, DE) ~ DP + TBL, data = kms, horizon = 1)`\
`#> `\
`#> Coefficients:`\
`#>          Estimate Std. Error t value Wald Ind Pr(> chi)    `\
`#> Ret:DP   0.006145   0.004557   1.349    1.819     0.177    `\
`#> DE:DP    0.311083   0.019396  16.039  257.244    <2e-16 ***`\
`#> Ret:TBL -0.080717   0.057701  -1.399    1.957     0.162    `\
`#> DE:TBL  -4.014921   0.277468 -14.470  209.377    <2e-16 ***`\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> Equation Wald statistics on 2 DF:`\
`#>      Wald p-value`\
`#> Ret 3.644  0.1617`\
`#> DE    488  <2e-16`\
`#> `\
`#> Joint Wald statistic:  493.9 on 4 DF, p-value < 2.2e-16`

Long horizons work as for
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md):

\
`s12`` ``<-`` `[`ivx_sys`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(``Ret``, ``DE``)`` ``~`` ``DP`` ``+`` ``TBL``, data ``=`` ``kms``, horizon ``=`` ``12``)`\
`s12``$``Wald_Eq`\
`#>        Ret         DE `\
`#>   3.998219 423.475294`

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the \\mr \times
mr\\ matrix of \\\mathrm{vec}(A)\\, column-major (all responses for the
first predictor, then the second, …), with names `response:predictor`;
custom restrictions can be tested from it directly.

\
`V`` ``<-`` `[`vcov`](https://rdrr.io/r/stats/vcov.html)`(``s``)`\
`a`` ``<-`` `[`as.vector`](https://rdrr.io/r/base/vector.html)`(`[`coef`](https://rdrr.io/r/stats/coef.html)`(``s``)``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``a``)`` ``<-`` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``V``)`\
`# is the effect of DP the same in both equations?`\
`H`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``-``1``, ``0``, ``0``)``, ``1``)`\
[`drop`](https://rdrr.io/r/base/drop.html)`(`[`t`](https://rdrr.io/r/base/t.html)`(``H`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``a``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` `[`solve`](https://rdrr.io/r/base/solve.html)`(``H`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``V`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` `[`t`](https://rdrr.io/r/base/t.html)`(``H``)``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``(``H`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``a``)``)`\
`#> [1] 233.1596`

## Validation

- For one response,
  [`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
  reproduces
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)
  exactly (tested for horizons 1 and 4), and the first equation’s Wald
  statistic in a system equals the univariate joint Wald statistic.
- Monte Carlo with two responses, a unit-root predictor and endogeneity
  correlations \\-0.95\\ and \\0.5\\, \\n = 500\\: empirical size of the
  joint 5% test 4.4%, of the equation tests 5.6% and 5.0%.

## Caveats

- The example uses `DE` as a second response only to illustrate the
  syntax; the theory assumes stationary responses.
- All equations share the same predictors and the same instrument;
  equation specific regressors are not supported.
- No Eicker–White option, no bootstrap and no `ivx_ar`/`ivx_ra`
  counterpart for systems.

## References

- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2023). Taking
  stock of long-horizon predictability tests: Are factor returns
  predictable? *Journal of Econometrics*, 237(2), 105380.
- Magdalinos, T. (2022). Least squares and IVX limit theory in systems
  of predictive regressions with GARCH innovations. *Econometric
  Theory*, 38(5), 875–912.

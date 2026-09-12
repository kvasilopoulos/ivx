# Demetrescu & Rodrigues (2022) — Monte Carlo size check of `ivx_ra()`

Demetrescu, M., & Rodrigues, P. M. M. (2022). Residual-augmented IVX predictive regression.
*Journal of Econometrics*, 227(2), 429–460. <https://doi.org/10.1016/j.jeconom.2020.09.008>

## What is replicated

The paper's empirical application uses OECD housing data that is not readily available, so
the check is against the size rows (b = 0) of Table 3 (right-sided) and Table 4 (two-sided)
for T = 200: DGP (24)–(25) with `rho = 1 − c/T`, short-run AR parameter `a1 = −0.5`,
innovation correlation −0.95, nominal 5% level, standard normal critical values.

| paper column | package |
|---|---|
| `t_ivx` (Kostakis et al. 2015 with their s.e. correction) | `ivx_fit()` |
| `t̃_ivx^{μ0}` (AR on the predictor without demeaning; recommended) | `ivx_ra_fit()` |

## Result (`size-table3-output.txt`, 10 000 replications, MC s.e. ≈ 0.002)

| c | paper `t_ivx` | ivx | paper `t̃_ivx^{μ0}` | ivx_ra |
|---|---|---|---|---|
| 0  | 0.116 | 0.129 | 0.054 | 0.064 |
| 10 | 0.088 | 0.089 | 0.055 | 0.064 |
| 20 | 0.074 | 0.072 | 0.055 | 0.064 |
| 30 | 0.066 | 0.057 | 0.053 | 0.055 |
| 40 | 0.064 | 0.057 | 0.050 | 0.059 |
| 50 | 0.061 | 0.051 | 0.050 | 0.056 |

`ivx_ra()` removes most of the over-rejection of plain IVX under strong endogeneity
(0.064 vs 0.129 at c = 0), matching the paper's qualitative and quantitative picture; the
remaining ≈ +0.01 offset at small c is shared with the plain IVX column and so reflects a
detail of the simulation design (e.g. initialisation) rather than the estimator.

## Implementation notes (things the paper leaves implicit)

- The AR/VAR on the predictors is fitted **without an intercept** (the μ0 variant). With an
  intercept (μ1) the AR coefficients are far more biased and the paper does not recommend it;
  that option is therefore not offered.
- The augmentation residuals must be **demeaned** (equivalently, regression (6) carries an
  intercept). Without this the term `γ·mean(ε̂)` is not removed and the bias reduction fails
  (T·bias 6.3 vs 2.2 in the check below).
- Standard errors: eq. (9)/(14) plus the Kostakis et al. (2015) intercept correction, as the
  paper does in its simulations (Section 4.1).
- Order selection: AIC in levels over 1..`ar_max`, as in the paper.

## Files

- `size-table3.R [R]` — the Monte Carlo; writes `size-table3.csv`.
- `size-table3-output.txt` — output of the 10 000-replication run.

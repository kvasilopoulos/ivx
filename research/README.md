# IVX literature — candidate extensions for the `ivx` package

`pdf/` holds the published journal version of each paper where one exists (arXiv / working-paper
versions otherwise; see "version" column); `txt/` holds `pdftotext -layout` extractions of the
same files (same basename). Collected 2026-09-12.

## What the package has today

- `ivx()` — Kostakis, Magdalinos & Stamatogiannis (2015, RFS) IVX Wald tests, short and long
  horizon, with the KMS finite-sample intercept correction. Tuning parameters are hard-coded in
  `src/ivx_fit_cpp.cpp`: `Rz = 1 - 1/n^0.95` (beta = 0.95, c_z = 1), Newey-West bandwidth
  `floor(n^(1/3))`, homoskedastic `sigma^2 * Z'Z` covariance.
- `ivx_ar()` — Yang, Long, Peng & Cai (2020, JASA) IVX-AR.
- Serial-correlation diagnostics (`ac_test_*`), `texreg` support.

## Papers

| file (pdf/ & txt/) | version | reference | why it matters for the package |
|---|---|---|---|
| `demetrescu-georgiev-rodrigues-taylor-2023-extensions-ivx` | journal (OA) | Demetrescu, Georgiev, Rodrigues & Taylor (2023) *Extensions to IVX methods of inference for return predictability*, JoE 237(2) 105271. [SD](https://www.sciencedirect.com/science/article/pii/S0304407622000586) · [BdP WP](https://ideas.repec.org/p/ptu/wpaper/w202104.html) | Eicker-White IVX covariance (Remarks 8-9), residual wild bootstrap (RWB) and fixed-regressor wild bootstrap (FRWB) algorithms (§4), subsample/episodic IVX tests (§3.2). One-sided t-tests. |
| `demetrescu-rodrigues-2022-residual-augmented-ivx` | journal | Demetrescu & Rodrigues (2022) *Residual-augmented IVX predictive regression*, JoE 227(2) 429-460. [SD](https://www.sciencedirect.com/science/article/abs/pii/S030440762030395X) · [BdP WP](https://ideas.repec.org/p/ptu/wpaper/w201605.html) | RA-IVX: 3-step bias-corrected IVX (AR(p) residuals of x → partial y on them → IVX). Estimator eq. (11), HC standard error eq. (12). |
| `hosseinkouchack-demetrescu-2021-finite-sample-ivx` | journal | Hosseinkouchack & Demetrescu (2021) *Finite-sample size control of IVX-based tests in predictive regressions*, ET 37 769-793. [CUP](https://www.cambridge.org/core/journals/econometric-theory/article/finitesample-size-control-of-ivxbased-tests-in-predictive-regressions/757376931B598A5196703802BCF4E2B3) | Convergence rate depends on IVX tuning (beta, c_z); higher-order bias terms; motivates exposing tuning params and finite-sample corrections. |
| `arxiv-2401.01064-robust-inference-multiple-predictive-regressions` | arXiv (unpublished) | Liao, Li & Fan (2024) *Robust inference for multiple predictive regressions with an application on bond risk premia*, arXiv:2401.01064. [arXiv](https://arxiv.org/abs/2401.01064) | Improved IVX (Q_m): sample-split weighted IV to remove the demeaning effect, bias & variance-enlargement corrections, LM-style restricted residuals. Algorithm 1 (p. 24). |
| `lee-2016-ivx-qr` | journal | Lee (2016) *Predictive quantile regression with persistent covariates: IVX-QR approach*, JoE 192(1) 105-118. [SD](https://www.sciencedirect.com/science/article/abs/pii/S0304407615003000) · [UCR pdf](https://economics.ucr.edu/wp-content/uploads/2019/10/Predictive-Quantile-Regression-with-Persistent-Covariates-IVX-QR-Approach.pdf) | IVX-QR. Practical test (§3.3): ordinary QR of y on the demeaned instrument z̃, self-normalised chi-square (Prop. 3.1); needs a sparsity estimate f_u(0). |
| `demetrescu-georgiev-rodrigues-taylor-2022-episodic-predictability` | journal | Demetrescu, Georgiev, Rodrigues & Taylor (2022) *Testing for episodic predictability in stock returns*, JoE 227(1) 85-113. [SD](https://www.sciencedirect.com/science/article/pii/S0304407620300026) · [BdP WP](https://www.bportugal.pt/sites/default/files/anexos/papers/wp201906.pdf) | Sup/ave/exp statistics over subsample IVX t-ratios ("pockets of predictability"), FRWB critical values. |
| `arxiv-2307.15151-predictability-tests-parameter-instability` | arXiv (unpublished) | Katsouris (2023) *Predictability tests robust against parameter instability*, arXiv:2307.15151. [arXiv](https://arxiv.org/abs/2307.15151) | Sup-Wald IVX tests for joint predictability + structural break; bootstrap critical values. |
| `katsouris-2023-structural-break-quantile-ivx` | arXiv (unpublished) | Katsouris (2023) *Structural break detection in quantile predictive regression models with persistent covariates*, arXiv:2302.05193. [arXiv](https://arxiv.org/abs/2302.05193) | Break tests in IVX-QR. |
| `katsouris-2023-bootstrap-ivx` | arXiv (unpublished) | Katsouris (2023) *Bootstrapping nonstationary autoregressive processes with predictive regression models*, arXiv:2307.14463. [arXiv](https://arxiv.org/abs/2307.14463) | Asymptotic validity of the bootstrap IVX estimator. |
| `katsouris-2024-cvar-doubly-ivx-qr` | arXiv (unpublished) | Katsouris (2024) *Estimating conditional Value-at-Risk with nonstationary quantile predictive regression models*, arXiv:2311.08218. [arXiv](https://arxiv.org/abs/2311.08218) | Doubly-IVX-corrected QR with generated regressors. |
| `arxiv-2309.14160-unified-dynamic-quantile-predictive` | arXiv (unpublished) | Katsouris (2023) *Unified inference for dynamic quantile predictive regression*, arXiv:2309.14160. [arXiv](https://arxiv.org/abs/2309.14160) | Dynamic (lagged-y) quantile predictive regression, IVX vs. alternatives. |
| `yang-liu-peng-cai-2021-unified-dynamic-predictive` | journal | Yang, Liu, Peng & Cai (2021) *Unified tests for a dynamic predictive regression*, JBES 39(3) 684-699. [IDEAS](https://ideas.repec.org/a/taf/jnlbes/v39y2021i3p684-699.html) | Companion to IVX-AR: tests with a lagged dependent variable. Not IVX-based (weighted empirical likelihood) but a natural neighbour of `ivx_ar`. |
| `kostakis-magdalinos-stamatogiannis-2023-long-horizon` | journal (OA) | Kostakis, Magdalinos & Stamatogiannis (2023) *Taking stock of long-horizon predictability tests: Are factor returns predictable?*, JoE 237(2) 105380. [SD](https://www.sciencedirect.com/science/article/pii/S0304407623000052) | Refined long-horizon IVX Wald; check current `horizon > 1` code against it. |
| `phillips-lee-2013-long-horizon-ivx` | journal | Phillips & Lee (2013) *Predictive regression under various degrees of persistence and robust long-horizon regression*, JoE 177(2) 250-264. [SD](https://doi.org/10.1016/j.jeconom.2013.04.011) | Original long-horizon IVX. |
| `xu-2020-multiple-horizon-predictability` | journal | Xu (2020) *Testing for multiple-horizon predictability: direct regression based versus implication based*, RFS 33(9) 4403-4443. [OUP](https://academic.oup.com/rfs/article-abstract/33/9/4403/5620729) | Implication-based multi-horizon test; IVX variant. |
| `magdalinos-2022-systems-ivx-garch` | journal | Magdalinos (2022) *Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations*, ET 38 875-912. [CUP](https://doi.org/10.1017/S0266466621000086) | Systems (multivariate y) IVX. |
| `liu-long-peng-yang-2024-unified-predictive-qr` | journal | Liu, Long, Peng & Yang (2024) *A unified inference for predictive quantile regression*, JASA 119(546). [T&F](https://doi.org/10.1080/01621459.2023.2203354) | Non-IVX alternative to IVX-QR. |

## Suggested extensions (ranked by value / effort)

1. **Expose IVX tuning** — `beta = 0.95`, `cz = 1`, NW bandwidth as arguments of `ivx()`/`ivx_fit()`
   (Hosseinkouchack & Demetrescu 2021; Lee 2016 uses beta = 0.95, cz = 5 for QR).
2. **Eicker-White covariance** — `vcov = c("iid", "hc")`: replace `sigma^2 * Z'Z` by `sum z z' u^2`
   (DGRT 2023, Remarks 8-9). Also report IVX t-ratios so one-sided tests are possible.
3. **Wild bootstrap p-values** — `ivx_boot(model, B, type = c("rwb", "frwb"))` following DGRT 2023
   Algorithms 4-5. Reuses `ivx_fit_cpp` in a loop and `auto_ar()` for the AR(p+1) on x (RWB).
4. **RA-IVX** — `ivx_ra()` implementing Demetrescu & Rodrigues (2022) eqs (11)-(12).
5. **IVX-QR** — `ivx_qr(formula, tau)` following Lee (2016) §3.3 (QR of y on z̃ via `quantreg::rq`,
   self-normalised chi-square). `quantreg` in Suggests.
6. **Episodic / subsample tests** — `ivx_episodic()` sup/ave IVX t-stats over windows with FRWB
   p-values (DGRT 2022; DGRT 2023 §3.2; Katsouris 2023 sup-Wald). Builds on (3).
7. **Liao-Li-Fan (2024) improved IVX** — `ivx(..., correction = "llf")` per Algorithm 1.
8. **Long-horizon audit** — verify `horizon > 1` against KMS (2023) and Phillips & Lee (2013).
9. **Systems IVX** — lift `stop("multivariate model is not available")` (Magdalinos 2022).

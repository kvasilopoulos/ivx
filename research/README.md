# IVX literature — candidate extensions for the `ivx` package

`pdf/` holds the published journal version of each paper where one exists (arXiv / working-paper
versions otherwise; see "version" column); `txt/` holds `pdftotext -layout` extractions of the
same files (same basename). Collected 2026-09-12.

See `replication/` for per-article replication scripts, data and results.

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
| `kostakis-magdalinos-stamatogiannis-2023-long-horizon` | journal (OA) | Kostakis, Magdalinos & Stamatogiannis (2023) *Taking stock of long-horizon predictability tests: Are factor returns predictable?*, JoE 237(2) 105380. [SD](https://www.sciencedirect.com/science/article/pii/S0304407623000052) | Refined long-horizon IVX Wald. Audited: `ivx(horizon = K)` implements eqs (15)/(23) exactly (K-sum response/regressor, single-lag instrument in the signal, K-sum instrument in the covariance, KMS intercept correction). |
| `phillips-lee-2013-long-horizon-ivx` | journal | Phillips & Lee (2013) *Predictive regression under various degrees of persistence and robust long-horizon regression*, JoE 177(2) 250-264. [SD](https://doi.org/10.1016/j.jeconom.2013.04.011) | Original long-horizon IVX. |
| `xu-2020-multiple-horizon-predictability` | journal | Xu (2020) *Testing for multiple-horizon predictability: direct regression based versus implication based*, RFS 33(9) 4403-4443. [OUP](https://academic.oup.com/rfs/article-abstract/33/9/4403/5620729) | Implication-based multi-horizon test; IVX variant. |
| `magdalinos-2022-systems-ivx-garch` | journal | Magdalinos (2022) *Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations*, ET 38 875-912. [CUP](https://doi.org/10.1017/S0266466621000086) | Systems (multivariate y) IVX. Implemented as `ivx_sys()` using the KMS (2023) Kronecker-form covariance. |
| `liu-long-peng-yang-2024-unified-predictive-qr` | journal | Liu, Long, Peng & Yang (2024) *A unified inference for predictive quantile regression*, JASA 119(546). [T&F](https://doi.org/10.1080/01621459.2023.2203354) | Non-IVX alternative to IVX-QR. |
| `campbell-yogo-2006-efficient-tests` | journal | Campbell & Yogo (2006) *Efficient tests of stock return predictability*, JFE 81(1) 27-60. [SD](https://doi.org/10.1016/j.jfineco.2005.05.008) | Bonferroni Q-test; package README to-do. Non-IVX benchmark. |
| `jansson-moreira-2006-optimal-inference` | journal | Jansson & Moreira (2006) *Optimal inference in regression models with nearly integrated regressors*, Econometrica 74(3) 681-714. [Wiley](https://doi.org/10.1111/j.1468-0262.2006.00679.x) | Conditional likelihood approach; package README to-do. Non-IVX benchmark. |
| `cavanagh-elliott-stock-1995-nearly-integrated` | journal | Cavanagh, Elliott & Stock (1995) *Inference in models with nearly integrated regressors*, ET 11(5) 1131-1147. [CUP](https://doi.org/10.1017/S0266466600009981) | Bonferroni method; package README to-do. |
| `elliott-2011-control-function` | journal | Elliott (2011) *A control function approach for testing the usefulness of trending variables in predictive regressions and econometric models*, JoE 164(1) 79-91. [SD](https://doi.org/10.1016/j.jeconom.2011.02.014) | Control function approach; package README to-do (cited there as "Elliott 2001"). |
| `demetrescu-2014-local-power-ivx` | journal | Demetrescu (2014) *Enhancing the local power of IVX-based tests in predictive regressions*, Econ. Letters 124 269-273. [SD](https://doi.org/10.1016/j.econlet.2014.05.032) | Lagged-y augmentation of the IVX regression; one `ivx()` argument. |
| `breitung-demetrescu-2015-iv-variable-addition` | journal | Breitung & Demetrescu (2015) *Instrumental variable and variable addition based inference in predictive regressions*, JoE 187 358-375. [SD](https://doi.org/10.1016/j.jeconom.2013.10.018) | Alternative IV choices / combined IVs, variable-addition tests. |
| `demetrescu-rodrigues-taylor-2023-transformed-long-horizon` | journal | Demetrescu, Rodrigues & Taylor (2023) *Transformed regression-based long-horizon predictability tests*, JoE 237 105316. [SD](https://doi.org/10.1016/j.jeconom.2022.06.006) | Competitor to KMS long-horizon; goes with the `horizon > 1` audit. |
| `fan-lee-2019-ivx-qr-garch` | SSRN WP (journal: JoE 213) | Fan & Lee (2019) *Predictive quantile regressions under persistence and conditional heteroskedasticity*, JoE 213 261-280. [SSRN](https://doi.org/10.2139/ssrn.3016449) | IVX-QR with GARCH errors; extends `ivx_qr()`. |
| `lee-shi-gao-2022-lasso-predictive` | journal | Lee, Shi & Gao (2022) *On LASSO for predictive regression*, JoE 229 322-349. [SD](https://doi.org/10.1016/j.jeconom.2021.02.002) | Many-predictor selection with IVX. |
| `amihud-hurvich-wang-2009-multiple-predictor` | journal | Amihud, Hurvich & Wang (2009) *Multiple-predictor regressions: hypothesis testing*, RFS 22(1) 413-434. [OUP](https://doi.org/10.1093/rfs/hhn056) | Multivariate augmented-regression (ARM); pairs with Elliott (2011). |
| `liu-yang-cai-peng-2019-unified-test` | journal | Liu, Yang, Cai & Peng (2019) *A unified test for predictability of asset returns regardless of properties of predicting variables*, JoE 208 141-159. [SD](https://doi.org/10.1016/j.jeconom.2018.09.009) | WEL test; same machinery as `ivx_ar`. |
| `harvey-leybourne-taylor-2021-simple-tests` | journal | Harvey, Leybourne & Taylor (2021) *Simple tests for stock return predictability with good size and power properties*, JoE 224 198-214. [SD](https://doi.org/10.1016/j.jeconom.2021.01.004) | Cheap OLS-based tests. |
| `chen-deo-yi-2013-uniform-inference` | journal | Chen, Deo & Yi (2013) *Uniform inference in predictive regression models*, JBES 31(4) 525-533. [T&F](https://doi.org/10.1080/07350015.2013.818008) | Restricted-likelihood test; pairs with Jansson & Moreira (2006). |
| `georgiev-harvey-leybourne-taylor-2018-parameter-instability` | journal | Georgiev, Harvey, Leybourne & Taylor (2018) *Testing for parameter instability in predictive regression models*, JoE 204(1) 101-118. [SD](https://doi.org/10.1016/j.jeconom.2018.01.005) | Pairs with the episodic / sup-Wald items. |
| `stock-1991-largest-root-ci` | journal | Stock (1991) *Confidence intervals for the largest autoregressive root in U.S. macroeconomic time series*, JME 28(3) 435-459. [SD](https://doi.org/10.1016/0304-3932(91)90034-L) | Inversion of a unit-root statistic into a CI for c; the Bonferroni step of Campbell & Yogo. |
| `elliott-rothenberg-stock-1996-dfgls` | journal | Elliott, Rothenberg & Stock (1996) *Efficient tests for an autoregressive unit root*, Econometrica 64(4) 813-836. [JSTOR](https://doi.org/10.2307/2171846) | DF-GLS statistic used by Campbell & Yogo and Harvey et al. |
| `hjalmarsson-2011-long-horizon` | journal | Hjalmarsson (2011) *New methods for inference in long-horizon regressions*, JFQA 46(3) 815-839. [CUP](https://doi.org/10.1017/S0022109011000135) | Long-horizon Bonferroni test; possible `horizon` for `cy_test()`. |

## Candidates not yet collected

None — all listed papers are in `pdf/`.

## Suggested extensions (ranked by value / effort)

Items already implemented are removed from this list and recorded in `NEWS.md`.

Deferred from the collected papers:

6. **Chen, Deo & Yi (2013)** — quasi-restricted-likelihood ratio test (WLSRL of the
   bivariate VAR, sup-bound critical value at c = 0 depending on the estimated innovation
   correlation). Deferred: the critical values are not tabulated in the paper and must be
   simulated from the local-to-unity limit functional for each delta; Remark 3 shows the
   right tail is close to chi-square(1) anyway, so the gain over `hlt_test()` / `ivx()` is
   small for the effort.
7. **Lee, Shi & Gao (2022)** — TAlasso variable selection; out of the package's scope
   (inference, not selection).

8. **Liao-Li-Fan (2024) improved IVX** — `ivx(..., correction = "llf")` per Algorithm 1.
   Deferred: unpublished, and Algorithm 1 (sample-split weights, bias term with its own
   tuning, variance-enlargement correction, LM residuals) needs the authors' code or the
   typeset equations to implement faithfully; not worth shipping from the text extraction.

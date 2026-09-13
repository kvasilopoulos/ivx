# ivx 1.2.0

* New vignette "Rolling IVX tests for bubble detection" showing how to build the
  rolling-window IVX test of Pavlidis, Paya & Peel (2017) from `ivx()` (#2).
* `ivx()` and `ivx_fit()` gain `beta`, `cz` and `bandwidth` arguments that expose
  the IVX instrument tuning (previously hard-coded to the Kostakis et al. (2015)
  values) and the Newey-West bandwidth.
* `ivx()` gains `robust = TRUE` for Eicker-White (heteroskedasticity-robust) IVX
  standard errors (Demetrescu, Georgiev, Rodrigues & Taylor, 2023).
* `ivx()` and `ivx_fit()` gain `lag_y = TRUE`: the lag-augmented IVX regression
  of Demetrescu (2014), which adds the lagged dependent variable (instrumented
  by itself) to raise local power under strong persistence and endogeneity.
  The joint Wald statistic tests the predictors only.
* `summary()` coefficient tables now report `Std. Error` and `t value` next to the
  individual Wald statistics; the fitted object stores `se` and `tstat`.
  Breaking: the table gains two columns, so code indexing `coef(summary(x))` by
  position must use column names (`"Wald Ind"`, `"Pr(> chi)"`) instead.
* New `ivx_boot()` implementing the residual wild bootstrap and fixed regressor
  wild bootstrap of Demetrescu et al. (2023), returning bootstrap p-values for
  the joint and individual Wald statistics and one-sided t-tests. Supports
  `cores > 1` via the \pkg{parallel} package; the regressor recursion of the
  residual wild bootstrap runs in C++.
* New `ivx_ra()` / `ivx_ra_fit()`: the residual-augmented (bias-reduced) IVX
  estimator of Demetrescu & Rodrigues (2022) with its heteroskedasticity-robust
  standard errors; returns an `ivx` object so `summary()`, `vcov()` etc. apply.
  `ivx_ra()` gains `horizon`: for `horizon > 1` it is the transformed-regression
  long-horizon test of Demetrescu, Rodrigues & Taylor (2023), which handles the
  overlap of the long-horizon regression without HAC estimation.
* New `ivx_qr()` / `ivx_qr_fit()`: the IVX-QR quantile predictability test of
  Lee (2016, Proposition 3.2) via \pkg{quantreg} (in Suggests); returns the
  estimated QR endogeneity `rho_tau` for the paper's tuning rule.
* New `ivx_qr_boot()`: moving block bootstrap percentile intervals and
  p-values for IVX-QR (Fan & Lee, 2019), robust to conditional
  heteroskedasticity and to the sparsity estimate.
* New `ivx_episodic()`: subsample (rolling, forward and backward recursive)
  IVX tests for pockets of predictability with sup/inf functionals and wild
  bootstrap p-values (Demetrescu et al. 2022, 2023 Section 3.2).
* Documented that the long-horizon statistic (`horizon > 1`) is the modified
  IVX-Wald of Kostakis, Magdalinos & Stamatogiannis (2023), eqs (15)/(23),
  after auditing the implementation against the paper.
* New `ivx_sys()` / `ivx_sys_fit()`: systems of predictive regressions with a
  matrix response (`cbind(y1, y2) ~ x`), short and long horizon, with the
  Kronecker-form IVX-Wald covariance of Kostakis et al. (2023); reports joint,
  per-equation and individual Wald statistics.
* New vignettes for each methodology (`ivx`, `ivx-sys`, `ivx-ar`, `ivx-ra`,
  `ivx-qr`, `robust-inference`, `ivx-episodic`) with the underlying
  statistics, replication results and caveats; pkgdown site reorganised
  (Bootstrap 5, MathJax rendering, grouped reference and articles).
* New `ivx_iv()` / `ivx_iv_fit()`: the 2SLS predictability tests of Breitung &
  Demetrescu (2015) with fractional-difference, long-difference and sine
  instruments and their recommended combination (`IVcomb`), Eicker-White
  standard errors.
* New `arm()` / `arm_fit()`: the multipredictor augmented regression method of
  Amihud, Hurvich & Wang (2009) - reduced-bias OLS with Nicholls-Pope corrected
  VAR(1) residuals as control variables and the paper's covariance estimator;
  a non-IVX benchmark for stationary persistent predictors.
* New `hlt_test()`: the hybrid switching t-test of Harvey, Leybourne & Taylor
  (2021) - standard or quasi-GLS-demeaned t-ratio with the paper's conservative
  critical values under strong persistence, normal critical values under weak
  persistence (ADF/MBIC switch).
* New `el_test()`: the unified empirical likelihood test of Liu, Yang, Cai &
  Peng (2019) for the predictive regression augmented with the lagged
  difference of the predictor; chi-square profile EL ratios whatever the
  persistence of the predictor, no tuning parameters.
* New `cy_test()`: the Bonferroni Q-test of Campbell & Yogo (2006) (the feasible
  Cavanagh, Elliott & Stock 1995 approach): DF-GLS confidence interval for the
  largest root inverted from simulated local-to-unity quantiles, Table 2 levels,
  Q-estimates with the AR(p) correction of Appendix A.
* New `elliott_cf()`: the control-function predictive regression of Elliott
  (2011) with user-supplied orthogonalising covariates and their lags, Wald
  test with Eicker-White standard errors and the remaining innovation
  correlation as a diagnostic.
* Fixed: weighted fits ignored the `horizon` argument; zero weights made
  `ivx()` fail (the dropped observations now get `NA` residuals and fitted
  values, with the same coefficients as fitting on the kept rows).
* Fixed: `ac_test_bg()` dispatched to the Box-Pierce method; `ac_test_lb()` and
  `ac_test_bp()` reported wrong p-values for non-consecutive `lag` vectors;
  the Breusch-Godfrey result now has the same `ac_test_` class and `pval`
  attribute as the other tests. `case.names()` no longer returns an empty vector.
* All formula interfaces share one model-frame routine, so `- 1` in a formula
  warns consistently and a matrix response is rejected consistently.
* `extract()` (texreg) now reports IVX standard errors instead of Wald statistics
  in the `se` slot.
* `ivx_ar()` and `ivx_ar_fit()` accept the same `beta`, `cz`, `bandwidth` and
  `robust` arguments as `ivx()`.

# ivx 1.1.1

* Patch version to fix minor issues.

# ivx 1.1.0

* Added `ivx_ar` that implements Yang, B., Long, W., Peng, L., & Cai, Z. (2020) 
new instrumental variable based  Wald statistic which accounts for serial 
correlation and heteroscedasticity in the error terms of the linear predictive regression model.
* Added the Yang et al. (2020) dataset named `ylpc`.
* Renamed the `monthly` and `quarterly` dataset into `kms` and `kms_quarterly`
* Removed dependency on `tibble` and `magrittr`.
* Added `texreg` functionality that converts regression output to LaTeX or HTML tables.
Specifically added `extract` methods for `ivx` and `ivx_ar`, which coefficients and GOF measures 
from a statistical object. 

# ivx 1.0.0

* Added a `NEWS.md` file to track changes to the package.

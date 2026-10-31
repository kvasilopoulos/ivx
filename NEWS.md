# ivx 2.0.0

## Breaking changes

* `summary()` coefficient tables now report `Std. Error` and `t value` next to the
  individual Wald statistics; the fitted object stores `se` and `tstat`.
  The table gains two columns, so code indexing `coef(summary(x))` by
  position must use column names (`"Wald Ind"`, `"Pr(> chi)"`) instead.

## Deprecated

* `drop1()`, `add1()` and `step()` on `ivx` objects. They compare
  residual sums of squares of an IV estimator, which is not a valid selection
  criterion; use the IVX Wald tests in `summary()`.

## New estimators and tests

* New `ivx_ra()` / `ivx_ra_fit()`: the residual-augmented (bias-reduced) IVX
  estimator of Demetrescu & Rodrigues (2022) with its heteroskedasticity-robust
  standard errors; returns an `ivx` object so `summary()`, `vcov()` etc. apply.
  `ivx_ra()` gains `horizon`: for `horizon > 1` it is the transformed-regression
  long-horizon test of Demetrescu, Rodrigues & Taylor (2023), which handles the
  overlap of the long-horizon regression without HAC estimation.
* New `ivx_qr()` / `ivx_qr_fit()`: the IVX-QR quantile predictability test of
  Lee (2016, Proposition 3.2) via `quantreg` (in Suggests); returns the
  estimated QR endogeneity `rho_tau` for the paper's tuning rule.
* New `ivx_qr_boot()`: moving block bootstrap percentile intervals and
  p-values for IVX-QR (Fan & Lee, 2019), robust to conditional
  heteroskedasticity and to the sparsity estimate.
* New `ivx_episodic()`: subsample (rolling, forward and backward recursive)
  IVX tests for pockets of predictability with sup/inf functionals and wild
  bootstrap p-values (Demetrescu et al. 2022, 2023 Section 3.2).
* New `ivx_sys()` / `ivx_sys_fit()`: systems of predictive regressions with a
  matrix response (`cbind(y1, y2) ~ x`), short and long horizon, with the
  Kronecker-form IVX-Wald covariance of Kostakis et al. (2023); reports joint,
  per-equation and individual Wald statistics.
* New `ivx_boot()` implementing the residual wild bootstrap and fixed regressor
  wild bootstrap of Demetrescu et al. (2023), returning bootstrap p-values for
  the joint and individual Wald statistics and one-sided t-tests. Supports
  `cores > 1` via the `parallel` package; the regressor recursion of the
  residual wild bootstrap runs in C++.
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

## New arguments and methods

* `ivx()` and `ivx_fit()` gain `beta`, `cz` and `bandwidth` arguments that expose
  the IVX instrument tuning (previously hard-coded to the Kostakis et al. (2015)
  values) and the Newey-West bandwidth.
* `ivx()` gains `robust = TRUE` for Eicker-White (heteroskedasticity-robust) IVX
  standard errors (Demetrescu, Georgiev, Rodrigues & Taylor, 2023).
* `ivx()` and `ivx_fit()` gain `lag_y = TRUE`: the lag-augmented IVX regression
  of Demetrescu (2014), which adds the lagged dependent variable (instrumented
  by itself) to raise local power under strong persistence and endogeneity.
  The joint Wald statistic tests the predictors only.
* `ivx_ar()` and `ivx_ar_fit()` accept the same `beta`, `cz`, `bandwidth` and
  `robust` arguments as `ivx()`.
* New `predict()` method for `ivx()` fits: the one-step-ahead conditional mean
  for the predictors in `newdata`, so coefficients estimated on a training
  sample can be projected out of sample (e.g., a fundamental price path).
* New `nobs()` method for `ivx` objects (and the extensions inheriting the class).
* `extract()` (texreg) now reports IVX standard errors instead of Wald statistics
  in the `se` slot.

## Bug fixes and output changes

* Fixed: `ivx_ar(ar = "auto")` ignored `ar_ic` and always chose the AR order
  by AIC; `ar_ic = "bic"` (the default) now selects by BIC. `ar = "forecast"`
  stops with an informative error when `forecast` is not installed.
* Arguments aligned with the papers they implement: `ivx_iv()` `d` is now the
  memory of the fractional instrument `(1 - L)^(1 - d) x[t-1]` and the long
  difference is `x[t-1] - x[t-k]` with `k <= t - 1`, as in
  Breitung & Demetrescu (2015, Section 2.2; the default `d = 1/2` is
  unchanged); `ivx_episodic()` `window` is the shortest window for every
  scheme and defaults to the settings of Demetrescu et al. (2022): 1/3 rolling,
  1/4 recursive; the `ivx_boot()` / `ivx_episodic()` `ar_max` default is the
  range of Demetrescu et al. (2023), `floor(4 * (T / 100)^(1/4)) + 1`;
  `elliott_cf()` defaults to the paper's conventional Wald test
  (`robust = FALSE`).
* Fixed: weighted fits ignored the `horizon` argument; zero weights made
  `ivx()` fail (the dropped observations now get `NA` residuals and fitted
  values, with the same coefficients as fitting on the kept rows).
* Fixed: `ac_test_bg()` dispatched to the Box-Pierce method; `ac_test_lb()` and
  `ac_test_bp()` reported wrong p-values for non-consecutive `lag` vectors;
  the Breusch-Godfrey result now has the same `ac_test_` class and `pval`
  attribute as the other tests. `case.names()` no longer returns an empty vector.
* `print()` of `ac_test()` on a numeric vector no longer errors (it indexed a
  Breusch-Godfrey column that only the `ivx` method computes).
* `ac_test()` and `ac_test_*()` now stop with an informative error on `ivx_ar`,
  `ivx_ra`, `ivx_qr`, `ivx_iv` and `arm` fits, which do not store the OLS
  residuals the tests use; previously they failed with an obscure error.
* `summary()` of `ivx_ar()` printed `Auto ()` without the information
  criterion; it now prints e.g. `Auto (BIC)`, and `ar = "forecast"` fits are
  labelled `forecast (BIC)` instead of `Fixed`.
* `print()` of `hlt_test()` and `cy_test()` splits its output into short lines
  (at most 70 characters).
* All formula interfaces share one model-frame routine, so `- 1` in a formula
  warns consistently and a matrix response is rejected consistently.
* The `ylpc` help page cited the wrong title, issue and pages for its source;
  it now cites Yang, Long, Peng & Cai (2020), JASA 115(532), 1598-1619.

## Documentation

* New vignette "Choosing a test: a decision workflow" (`vignette("workflow")`)
  that walks from the plain `ivx()` fit through the diagnostics (`ac_test()`,
  `delta()`, ARCH check) to the extension each one calls for.
* New vignette "Rolling IVX tests for bubble detection" showing how to build the
  rolling-window IVX test of Pavlidis, Paya & Peel (2017) from `ivx()` (#2).
* New vignettes for each methodology (`ivx`, `ivx-sys`, `ivx-ar`, `ivx-ra`,
  `ivx-qr`, `robust-inference`, `ivx-episodic`) with the underlying
  statistics, replication results and caveats; pkgdown site reorganised
  (Bootstrap 5, MathJax rendering, grouped reference and articles).
* Documented that the long-horizon statistic (`horizon > 1`) is the modified
  IVX-Wald of Kostakis, Magdalinos & Stamatogiannis (2023), eqs (15)/(23),
  after auditing the implementation against the paper.

# ivx 1.1.1

* Fixed `ivx_ar()` when `auto.arima()` selects the AR order (#3).
* Added an S3 method for `ac_test_bg()`, and `extract()` for `ivx` now checks
  that `texreg` is installed.
* Removed the C++11 specification from `src/Makevars` and `src/Makevars.win`
  to address a CRAN NOTE; `inst/CITATION` now uses `bibentry()` instead of the
  old-style `citEntry()`.
* Documentation and pkgdown fixes.

# ivx 1.1.0

* Added `ivx_ar()`, which implements the instrumental variable based Wald
  statistic of Yang, B., Long, W., Peng, L., & Cai, Z. (2020). It accounts for
  serial correlation and heteroscedasticity in the error terms of the linear
  predictive regression model.
* Added the Yang et al. (2020) dataset named `ylpc`.
* Breaking: renamed the `monthly` and `quarterly` datasets to `kms` and
  `kms_quarterly`.
* Removed the dependency on `tibble` and `magrittr`.
* Added `texreg` support that converts regression output to LaTeX or HTML
  tables: `extract()` methods for `ivx` and `ivx_ar` extract coefficients and
  GOF measures from the fitted object.

# ivx 1.0.0

* Added a `NEWS.md` file to track changes to the package.

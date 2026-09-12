# Changelog

## ivx 1.2.0

- New vignette “Rolling IVX tests for bubble detection” showing how to
  build the rolling-window IVX test of Pavlidis, Paya & Peel (2017) from
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)
  ([\#2](https://github.com/kvasilopoulos/ivx/issues/2)).
- [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) and
  [`ivx_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_fit.md)
  gain `beta`, `cz` and `bandwidth` arguments that expose the IVX
  instrument tuning (previously hard-coded to the Kostakis et al. (2015)
  values) and the Newey-West bandwidth.
- [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) gains
  `robust = TRUE` for Eicker-White (heteroskedasticity-robust) IVX
  standard errors (Demetrescu, Georgiev, Rodrigues & Taylor, 2023).
- [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) and
  [`ivx_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_fit.md)
  gain `lag_y = TRUE`: the lag-augmented IVX regression of Demetrescu
  (2014), which adds the lagged dependent variable (instrumented by
  itself) to raise local power under strong persistence and endogeneity.
  The joint Wald statistic tests the predictors only.
- [`summary()`](https://rdrr.io/r/base/summary.html) coefficient tables
  now report `Std. Error` and `t value` next to the individual Wald
  statistics; the fitted object stores `se` and `tstat`. Breaking: the
  table gains two columns, so code indexing `coef(summary(x))` by
  position must use column names (`"Wald Ind"`, `"Pr(> chi)"`) instead.
- New
  [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  implementing the residual wild bootstrap and fixed regressor wild
  bootstrap of Demetrescu et al. (2023), returning bootstrap p-values
  for the joint and individual Wald statistics and one-sided t-tests.
  Supports `cores > 1` via the package; the regressor recursion of the
  residual wild bootstrap runs in C++.
- New
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  /
  [`ivx_ra_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra_fit.md):
  the residual-augmented (bias-reduced) IVX estimator of Demetrescu &
  Rodrigues (2022) with its heteroskedasticity-robust standard errors;
  returns an `ivx` object so
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) etc. apply.
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  gains `horizon`: for `horizon > 1` it is the transformed-regression
  long-horizon test of Demetrescu, Rodrigues & Taylor (2023), which
  handles the overlap of the long-horizon regression without HAC
  estimation.
- New
  [`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
  /
  [`ivx_qr_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_fit.md):
  the IVX-QR quantile predictability test of Lee (2016, Proposition 3.2)
  via (in Suggests); returns the estimated QR endogeneity `rho_tau` for
  the paper’s tuning rule.
- New
  [`ivx_qr_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_boot.md):
  moving block bootstrap percentile intervals and p-values for IVX-QR
  (Fan & Lee, 2019), robust to conditional heteroskedasticity and to the
  sparsity estimate.
- New
  [`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md):
  subsample (rolling, forward and backward recursive) IVX tests for
  pockets of predictability with sup/inf functionals and wild bootstrap
  p-values (Demetrescu et al. 2022, 2023 Section 3.2).
- Documented that the long-horizon statistic (`horizon > 1`) is the
  modified IVX-Wald of Kostakis, Magdalinos & Stamatogiannis (2023), eqs
  (15)/(23), after auditing the implementation against the paper.
- New
  [`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
  /
  [`ivx_sys_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys_fit.md):
  systems of predictive regressions with a matrix response
  (`cbind(y1, y2) ~ x`), short and long horizon, with the Kronecker-form
  IVX-Wald covariance of Kostakis et al. (2023); reports joint,
  per-equation and individual Wald statistics.
- New vignettes for each methodology (`ivx`, `ivx-sys`, `ivx-ar`,
  `ivx-ra`, `ivx-qr`, `robust-inference`, `ivx-episodic`) with the
  underlying statistics, replication results and caveats; pkgdown site
  reorganised (Bootstrap 5, MathJax rendering, grouped reference and
  articles).
- New
  [`ivx_iv()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md)
  /
  [`ivx_iv_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv_fit.md):
  the 2SLS predictability tests of Breitung & Demetrescu (2015) with
  fractional-difference, long-difference and sine instruments and their
  recommended combination (`IVcomb`), Eicker-White standard errors.
- New [`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md) /
  [`arm_fit()`](https://kvasilopoulos.github.io/ivx/reference/arm_fit.md):
  the multipredictor augmented regression method of Amihud, Hurvich &
  Wang (2009) - reduced-bias OLS with Nicholls-Pope corrected VAR(1)
  residuals as control variables and the paper’s covariance estimator; a
  non-IVX benchmark for stationary persistent predictors.
- Fixed: weighted fits ignored the `horizon` argument.
- [`extract()`](https://magrittr.tidyverse.org/reference/aliases.html)
  (texreg) now reports IVX standard errors instead of Wald statistics in
  the `se` slot.
- [`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
  and
  [`ivx_ar_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar_fit.md)
  accept the same `beta`, `cz`, `bandwidth` and `robust` arguments as
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md).

## ivx 1.1.1

CRAN release: 2025-09-20

- Patch version to fix minor issues.

## ivx 1.1.0

CRAN release: 2020-11-24

- Added `ivx_ar` that implements Yang, B., Long, W., Peng, L., &
  Cai, Z. (2020) new instrumental variable based Wald statistic which
  accounts for serial correlation and heteroscedasticity in the error
  terms of the linear predictive regression model.
- Added the Yang et al. (2020) dataset named `ylpc`.
- Renamed the `monthly` and `quarterly` dataset into `kms` and
  `kms_quarterly`
- Removed dependency on `tibble` and `magrittr`.
- Added `texreg` functionality that converts regression output to LaTeX
  or HTML tables. Specifically added `extract` methods for `ivx` and
  `ivx_ar`, which coefficients and GOF measures from a statistical
  object.

## ivx 1.0.0

CRAN release: 2019-05-04

- Added a `NEWS.md` file to track changes to the package.

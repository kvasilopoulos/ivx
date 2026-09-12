# ivx 1.2.0

* New vignette "Rolling IVX tests for bubble detection" showing how to build the
  rolling-window IVX test of Pavlidis, Paya & Peel (2017) from `ivx()` (#2).
* `ivx()` and `ivx_fit()` gain `beta`, `cz` and `bandwidth` arguments that expose
  the IVX instrument tuning (previously hard-coded to the Kostakis et al. (2015)
  values) and the Newey-West bandwidth.
* `ivx()` gains `robust = TRUE` for Eicker-White (heteroskedasticity-robust) IVX
  standard errors (Demetrescu, Georgiev, Rodrigues & Taylor, 2023).
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
* Fixed: weighted fits ignored the `horizon` argument.
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

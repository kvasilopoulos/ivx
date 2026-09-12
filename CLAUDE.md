# ivx

R package: IVX (extended instrumental variable) inference for predictive
regressions with persistent regressors. Rcpp/Armadillo backend.

## Map

- `R/ivx.R` —
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) /
  [`ivx_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_fit.md)
  (Kostakis, Magdalinos & Stamatogiannis 2015); core is
  `src/ivx_fit_cpp.cpp`. Tuning via `beta`, `cz`, `bandwidth`;
  `robust = TRUE` for Eicker-White; `lag_y = TRUE` lag-augmented
  (Demetrescu 2014).
- `R/ivx-ar.R` —
  [`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
  IVX-AR (Yang, Long, Peng & Cai 2020).
- `R/ivx-ra.R` —
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  residual-augmented IVX (Demetrescu & Rodrigues 2022); `horizon > 1` is
  the transformed-regression test (Demetrescu, Rodrigues & Taylor 2023).
- `R/ivx-qr.R` —
  [`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
  IVX quantile regression (Lee 2016),
  [`ivx_qr_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_boot.md)
  block bootstrap (Fan & Lee 2019); uses `quantreg`.
- `R/ivx-sys.R` —
  [`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
  systems IVX (Magdalinos 2022); `src/ivx_sys_fit_cpp.cpp`.
- `R/ivx-boot.R` —
  [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  residual / fixed-regressor wild bootstrap (DGRT 2023);
  `src/var_sim.cpp` for the regressor recursion.
- `R/ivx-episodic.R` —
  [`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md)
  subsample sup/ave tests (DGRT 2022); `src/sub_ivx.cpp`.
- `R/ivx-iv.R` —
  [`ivx_iv()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md)
  2SLS with sine/fractional/long-difference instruments (Breitung &
  Demetrescu 2015).
- Non-IVX benchmarks: `R/arm.R`
  ([`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md),
  Amihud-Hurvich-Wang 2009), `R/hlt-test.R`
  ([`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md),
  Harvey-Leybourne-Taylor 2021), `R/el-test.R`
  ([`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md),
  Liu-Yang-Cai-Peng 2019; `el_ratio()` is a generic Owen EL solver).
- `R/ac_test.R` — serial-correlation diagnostics. `R/methods*.R`,
  `R/extract-texreg-methods.R` — S3 / texreg support. `R/auto-ar.R` — AR
  order selection helper.
- `tests/testthat/test-<file>.R` mirrors each `R/` file.
- `research/` — literature for extensions (not built). See
  `research/CLAUDE.md`.

## Git

- Never add a `Co-Authored-By` trailer or any AI-attribution line to
  commit messages.

## Changelog

- Every user-visible change (new function/argument, changed output, bug
  fix, breaking change) gets a bullet in `NEWS.md` under the current
  development version, in the same commit as the change.

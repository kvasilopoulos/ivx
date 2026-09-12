# ivx

R package: IVX (extended instrumental variable) inference for predictive
regressions with persistent regressors. Rcpp/Armadillo backend.

## Map

- `R/ivx.R` —
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) /
  [`ivx_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_fit.md)
  (Kostakis, Magdalinos & Stamatogiannis 2015); core is
  `src/ivx_fit_cpp.cpp`. Tuning via `beta`, `cz`, `bandwidth`;
  `robust = TRUE` for Eicker-White.
- `R/ivx-ar.R` —
  [`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
  IVX-AR (Yang, Long, Peng & Cai 2020).
- `R/ivx-ra.R` —
  [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  residual-augmented IVX (Demetrescu & Rodrigues 2022).
- `R/ivx-qr.R` —
  [`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
  IVX quantile regression (Lee 2016); uses `quantreg`.
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

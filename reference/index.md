# Package index

## IVX estimation

Predictive regressions with persistent regressors, short and long
horizon (Kostakis, Magdalinos & Stamatogiannis 2015, 2023).

- [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)
  [`print(`*`<ivx>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)
  : Fitting IVX Models
- [`ivx_sys()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
  [`print(`*`<ivx_sys>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
  [`summary(`*`<ivx_sys>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys.md)
  : Fitting Systems of IVX Predictive Regressions

## Extensions

Modifications of the IVX estimator for specific error structures or
targets.

- [`ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
  [`print(`*`<ivx_ar>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)
  : Fitting IVX-AR Models
- [`ivx_ra()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  [`print(`*`<ivx_ra>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra.md)
  : Fitting Residual-Augmented IVX Models
- [`ivx_qr()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
  [`print(`*`<ivx_qr>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
  [`summary(`*`<ivx_qr>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr.md)
  : Fitting IVX Quantile Predictive Regressions
- [`ivx_iv()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md)
  [`print(`*`<ivx_iv>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv.md)
  : IV Predictability Tests with Combined Instruments

## Robust and subsample inference

Wild bootstrap p-values and tests for pockets of predictability
(Demetrescu et al. 2022, 2023).

- [`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  [`print(`*`<ivx_boot>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
  : Wild Bootstrap Inference for IVX Models
- [`ivx_qr_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_boot.md)
  [`print(`*`<ivx_qr_boot>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_boot.md)
  : Moving Block Bootstrap for IVX-QR
- [`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md)
  [`print(`*`<ivx_episodic>`*`)`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md)
  : Subsample IVX Tests for Episodic Predictability

## Non-IVX benchmarks

Reduced-bias OLS alternatives for stationary persistent predictors.

- [`arm()`](https://kvasilopoulos.github.io/ivx/reference/arm.md)
  [`print(`*`<arm>`*`)`](https://kvasilopoulos.github.io/ivx/reference/arm.md)
  : Augmented Regression Method (Amihud, Hurvich & Wang)
- [`hlt_test()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md)
  [`hlt_test_fit()`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md)
  [`print(`*`<hlt_test>`*`)`](https://kvasilopoulos.github.io/ivx/reference/hlt_test.md)
  : Hybrid t-Test for Return Predictability (Harvey, Leybourne & Taylor)
- [`el_test()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md)
  [`el_test_fit()`](https://kvasilopoulos.github.io/ivx/reference/el_test.md)
  [`print(`*`<el_test>`*`)`](https://kvasilopoulos.github.io/ivx/reference/el_test.md)
  : Unified Empirical Likelihood Test for Predictability (Liu, Yang, Cai
  & Peng)

## Fitter functions

Low-level functions taking a response and a design matrix.

- [`ivx_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_fit.md)
  [`ivx_wfit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_fit.md)
  : Fitter Functions for IVX Models
- [`ivx_ar_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar_fit.md)
  : Fitter Functions for IVX-AR Models
- [`ivx_ra_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_ra_fit.md)
  : Fitter Function for Residual-Augmented IVX Models
- [`ivx_qr_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_qr_fit.md)
  : Fitter Function for IVX-QR Models
- [`arm_fit()`](https://kvasilopoulos.github.io/ivx/reference/arm_fit.md)
  : Fitter Function for the Augmented Regression Method
- [`ivx_iv_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_iv_fit.md)
  : Fitter Function for IV Predictability Tests
- [`ivx_sys_fit()`](https://kvasilopoulos.github.io/ivx/reference/ivx_sys_fit.md)
  : Fitter Function for Systems IVX Models

## Serial correlation tests

- [`ac_test()`](https://kvasilopoulos.github.io/ivx/reference/ac_test.md)
  : Autocorrelation tests
- [`ac_test_wald()`](https://kvasilopoulos.github.io/ivx/reference/ac_test_.md)
  [`ac_test_lb()`](https://kvasilopoulos.github.io/ivx/reference/ac_test_.md)
  [`ac_test_bp()`](https://kvasilopoulos.github.io/ivx/reference/ac_test_.md)
  [`ac_test_bg()`](https://kvasilopoulos.github.io/ivx/reference/ac_test_.md)
  : Tests for autocorrelation

## Methods

- [`summary(`*`<ivx>`*`)`](https://kvasilopoulos.github.io/ivx/reference/summary.ivx.md)
  [`print(`*`<summary.ivx>`*`)`](https://kvasilopoulos.github.io/ivx/reference/summary.ivx.md)
  : Summarizing IVX Model Fits

- [`summary(`*`<ivx_ar>`*`)`](https://kvasilopoulos.github.io/ivx/reference/summary.ivx_ar.md)
  [`print(`*`<summary.ivx_ar>`*`)`](https://kvasilopoulos.github.io/ivx/reference/summary.ivx_ar.md)
  : Summarizing IVX-AR Model Fits

- [`vcov(`*`<ivx>`*`)`](https://kvasilopoulos.github.io/ivx/reference/vcov.ivx.md)
  [`vcov(`*`<summary.ivx>`*`)`](https://kvasilopoulos.github.io/ivx/reference/vcov.ivx.md)
  : Calculate Variance-Covariance Matrix for a Fitted Model Object

- [`delta()`](https://kvasilopoulos.github.io/ivx/reference/delta.md) :
  Calculate the delta coefficient

- [`extract.ivx()`](https://kvasilopoulos.github.io/ivx/reference/extract.ivx.md)
  [`extract.ivx_ar()`](https://kvasilopoulos.github.io/ivx/reference/extract.ivx.md)
  :

  `extract` method for `ivx` objects

## Data

- [`kms`](https://kvasilopoulos.github.io/ivx/reference/kms.md) : KMS
  Monthly data
- [`kms_quarterly`](https://kvasilopoulos.github.io/ivx/reference/kms_quarterly.md)
  : KMS Quarterly data
- [`monthly`](https://kvasilopoulos.github.io/ivx/reference/monthly.md)
  : Monthly dataset of KMS
- [`quarterly`](https://kvasilopoulos.github.io/ivx/reference/quarterly.md)
  : Quarterly dataset of KMS
- [`ylpc`](https://kvasilopoulos.github.io/ivx/reference/ylpc.md) : YLPC
  Quarterly data

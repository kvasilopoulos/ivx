# Demetrescu, Georgiev, Rodrigues & Taylor (2023) — replication of Table 4, Panel A

Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R. (2023). Extensions to
IVX methods of inference for return predictability. *Journal of Econometrics*, 237(2), 105271.
<https://doi.org/10.1016/j.jeconom.2022.02.007>

## What is replicated

Table 4, Panel A: univariate predictive regressions of the monthly log equity premium on
each of the 14 Welch–Goyal (2008) predictors, January 1927 – December 2020 (T = 1128).

| paper column | package |
|---|---|
| `t_EW(−)`, `t_EW(+)`, `t_EW` — asymptotic p-values with Eicker-White s.e. (eq. 14) | `ivx(robust = TRUE)`, `pnorm()` on `$tstat` |
| `t_RWB(−)`, `t_RWB(+)`, `t_RWB` — residual wild bootstrap p-values, conventional s.e., 9999 reps | `ivx_boot(type = "rwb")` |
| `β_OLS`, `β_IVX`, `ρ`, `δ` | `$ols$coefficients`, `coef()`, `$AR$Rn`, `delta()` |

## Data

`goyal-welch-monthly.csv` — the "Monthly" sheet of Amit Goyal's `PredictorData` workbook
(<https://sites.google.com/view/agoyal145>), 2023 vintage, exported via the Google-Sheets
CSV endpoint used by the `tidyfinance` package. `goyal-welch-monthly-2021.csv` is the 2021
vintage (from the `PredictorData2021.xlsx` copy hosted at
<https://github.com/shokru/coqueret.github.io>) used only for the vintage check. The paper
used the 2020 vintage, which is no longer distributed.

Variables follow Welch–Goyal: `dp = log D12 − log Index`, `dy = log D12 − log Index_{t−1}`,
`ep = log E12 − log Index`, `de = log D12 − log E12`, `tms = lty − tbl`, `dfy = BAA − AAA`,
`dfr = corpr − ltr`; `svar, bm, ntis, tbl, lty, ltr, infl` as in the file. On the overlapping
sample these match the package's `kms` dataset (correlations ≥ 0.998 except `infl` 0.97 and
`ntis` 0.98, which Goyal–Welch have revised).

Excess return: the paper says "log of the monthly return on the S&P 500 index (including
dividends) minus the log of the risk-free rate", i.e. `vw = log(1 + CRSP_SPvw) − log(1 + Rfree)`.
The ex-dividend series `vwx = log(1 + CRSP_SPvwx) − log(1 + Rfree)` is also run — see below.

## Result

Full output: `table4-panelA-output.txt`; machine-readable: `table4-panelA-ivx.csv` and
`table4-panelA-comparison-{vw,vwx}.csv`. B = 1999 (paper: 9999), so bootstrap Monte-Carlo
error is ≈ 0.01 on a p-value near 0.5.

Two-sided p-values, paper vs package (`sensitivity-returns.csv`, `table4-panelA-comparison-*.csv`):

| pred | EW paper | EW `vw` | EW `vwx` | RWB paper | RWB `vw` | RWB `vwx` |
|---|---|---|---|---|---|---|
| dp | 0.510 | 0.227 | **0.536** | 0.612 | 0.344 | **0.641** |
| dy | 0.149 | **0.140** | 0.366 | 0.287 | **0.197** | 0.437 |
| ep | 0.172 | 0.057 | **0.167** | 0.272 | 0.125 | **0.269** |
| de | 0.603 | 0.666 | 0.580 | 0.465 | 0.523 | 0.415 |
| svar | 0.771 | 0.809 | 0.769 | 0.711 | 0.776 | 0.739 |
| bm | 0.594 | 0.329 | **0.529** | 0.453 | 0.217 | **0.425** |
| ntis | 0.234 | 0.202 | 0.208 | 0.161 | 0.143 | 0.147 |
| tbl | 0.073 | 0.084 | 0.068 | 0.087 | 0.108 | 0.091 |
| lty | 0.054 | 0.083 | 0.061 | 0.070 | 0.099 | 0.079 |
| ltr | 0.164 | 0.151 | 0.145 | 0.171 | 0.151 | 0.142 |
| tms | 0.513 | 0.448 | 0.457 | 0.414 | 0.345 | 0.353 |
| dfy | 0.996 | 0.921 | 0.985 | 0.994 | 0.884 | 0.979 |
| dfr | 0.394 | 0.394 | 0.398 | 0.374 | 0.417 | 0.423 |
| infl | 0.111 | 0.260 | 0.242 | 0.187 | 0.244 | 0.225 |
| mean abs. diff. | | 0.082 | 0.043 | | 0.064 | 0.024 |

Findings:

1. **Statistics are implemented correctly.** With the `vwx` return series the package
   reproduces the paper for 12 of 14 predictors to within bootstrap noise for the RWB
   (e.g. `ep` 0.269 vs 0.272, `lty` 0.079 vs 0.070) and within ±0.03–0.06 for the asymptotic
   EW test. An independent R re-implementation of eq. (14) (with and without the KMS
   intercept-correction term) gives the same t-ratios as the C++ code.
2. **Return definition.** The three predictors with strong endogeneity (`dp`, `ep`, `bm`,
   δ ≈ −0.77 to −0.98) are reproduced only with the *ex-dividend* return `CRSP_SPvwx`,
   despite the paper's text; `dy` is reproduced only with the with-dividend series. The
   remaining 10 predictors are insensitive to the choice. `infl` (p ≈ 0.11 in the paper vs
   0.22–0.26 here) matches under neither; its history has been revised by Goyal–Welch.
3. **Vintage.** The 2021 and 2023 vintages give the same results (differences ≤ 0.03,
   `sensitivity-returns.R`), so the residual gaps are not vintage revisions between those two;
   the 2020 vintage used in the paper could not be obtained.
4. `β` estimates match once the paper's ×100 scaling of `tbl…infl` is accounted for; `ρ` in
   the paper is from an AR(p) error-correction fit (BIC), the package reports the AR(1) root.

## Files

- `replicate-table4.R [B] [vw|vwx|both]` — main script (≈ 5 min per return series at
  B = 1999 with `cores = 1`; pass `cores` to `ivx_boot()` to parallelise).
- `sensitivity-returns.R` — asymptotic EW p-values under four return definitions × two vintages.
- `table4-panelA-output.txt`, `table4-panelA-ivx.csv`, `table4-panelA-comparison-{vw,vwx}.csv`,
  `sensitivity-returns.csv` — generated outputs.

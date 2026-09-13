# Monte Carlo checks for the 1.2.0 extensions

`size-checks.R` reproduces, from one script with fixed seeds, the rejection frequencies
quoted in the vignettes for the methods added in version 1.2.0. Default 500 replications
per cell (`Rscript size-checks.R 2000` for tighter numbers); output in
`size-checks-output.txt`. Monte Carlo s.e. at 5% with 500 replications is about 0.01.

| block | function | design | paper comparison | result |
|---|---|---|---|---|
| 1 | `ivx_ra(horizon = h)` | DRT (2023): T = 500, δ = −0.95, c ∈ {0, −20}, h ∈ {10, 20}, two-sided 5% | paper range [0.023, 0.058] | 0.028–0.056 |
| 2 | `ivx(lag_y = TRUE)` | Demetrescu (2014) Table 1: T = 100, δ = −0.98, η = 0, β = b/T | b = 0: 4.6 / 7.9; b = 5: 8.3 / 63.9 (IVX / augmented) | 5.6 / 4.2; 5.6 / 58.6 |
| 3 | `ivx_qr_boot()` | Fan & Lee (2019) study 2: ARCH(1) α₁ = 0.9, ρ = −0.9, c = 0, n = 200, τ = 0.1 | asymptotic test oversized in the tails, MBB near nominal | 20.8% vs 7.8% |
| 4 | `arm()` | AHW (2009) Case 1, n = 100 | 6–10% (n = 50), 5–7% (n = 200) | t: 5.8 / 8.0, Wald 7.4 |
| 5 | `ivx_iv()` | B&D (2015) Table 1: T = 250, corr 0.9, ρ = 1, 10% level, b ∈ {0, 10, 20} | comb 11.2 / 65.7 / 91.3; sin 9.9 / 61.4 / 79.5; frac 11.1 / 33.4 / 66.9; diff 12.5 / 33.6 / 66.1 | comb 13.2 / 52.4 / 85.4; sin 10.6 / 38.8 / 56.8; frac 12.0 / 44.8 / 78.4; diff 11.4 / 34.2 / 66.2 |
| 6 | `hlt_test()` | HLT (2021): T = 200, upper-tail 5%, ρ_xy ∈ {−0.9, 0.9}, ρ ∈ {1, 0.95, 0.5} | "very little deviation from nominal size, apart from some undersize for positive ρ_xy in the more persistent cases" | 7.6 / 5.0 / 4.4 and 0.2 / 2.8 / 6.0 |
| 7 | `el_test()` | LYCP (2019), Gaussian innovations: n = 500, β₁ = 0.5, β₂ = 0, ρ ∈ {1, 0.9} | "accurate size for H₀: β₂ = 0" | 7.6 / 5.2 |
| 8 | `cy_test()` | CY (2006): T = 200, δ ∈ {−0.9, −0.5}, c ∈ {0, −20}, one-sided 5% | right-tailed ≈ 4–5%, left-tailed as low as 1.2% (Section 3.4) | right 4.0–5.0, left 1.0–2.4 |
| 9 | `elliott_cf()` | Elliott (2011) idea: z = v + noise, δ = −0.9, c = 0, T = 200 | covariate restores nominal size | 4.0% vs 23.4% for OLS |

Sizes and orderings match the papers. The one systematic gap is block 5: the sine
instrument has less power and the fractional instrument more than in the paper's table;
the paper does not report how the deterministic instrument is demeaned or initialised
(see the `ivx-iv` vignette).

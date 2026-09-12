# Kostakis, Magdalinos & Stamatogiannis (2015) — replication

Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). Robust econometric inference
for stock return predictability. *Review of Financial Studies*, 28(5), 1506–1553.
<https://doi.org/10.1093/rfs/hhu139>

## What is replicated

Package function: `ivx()` with default tuning (`beta = 0.95`, `cz = 1`, Newey-West
bandwidth `n^(1/3)`, homoskedastic covariance with the KMS intercept correction).

| paper | regression | quantities |
|---|---|---|
| Table 6 (p. 1531) | `Ret ~ DE` | IVX slope, individual Wald, δ |
| Table 8 (p. 1537) | `Ret ~ DP + TBL` | IVX slopes, joint Wald |
| Table 11 (p. 1544) | `Ret ~ DE`, horizons 4–60 | individual Wald |
| Table 13 (p. 1547) | `Ret ~ EP + TBL`, horizons 4–60 | individual and joint Wald |

Data: `kms` dataset in the package (the authors' monthly file, 1926-12 to 2012-12).

## Result

All published numbers are reproduced to the printed precision (3–4 decimals);
see the script output. The same assertions live in `tests/testthat/test-ivx.R`,
so they are checked on every test run.

## Files

- `replicate-tables.R` — prints paper vs package values side by side.

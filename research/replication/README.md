# Replications

One folder per article. Each folder has a `README.md` (what is replicated, data source,
outcome), the script(s), any data not shipped with the package, and the generated
output files. Scripts resolve paths relative to their own location and load the package
from the repository root with `pkgload::load_all()`, so they can be run from anywhere:

```sh
Rscript research/replication/<article>/<script>.R
```

| folder | article | package feature | status |
|---|---|---|---|
| `kostakis-magdalinos-stamatogiannis-2015/` | KMS (2015, RFS) Tables 6, 8, 11, 13 | `ivx()` short and long horizon | exact to printed precision |
| `demetrescu-georgiev-rodrigues-taylor-2023/` | DGRT (2023, JoE) Table 4 Panel A | `ivx(robust = TRUE)`, `ivx_boot(type = "rwb")` | 12/14 predictors within bootstrap noise (ex-dividend returns); see folder README |
| `demetrescu-rodrigues-2022/` | DR (2022, JoE) Table 3 size (T = 200, b = 0) | `ivx_ra()` | within MC error; see folder README |
| `extensions-2026/` | size/power checks for every method added in 1.2.0 (DRT 2023, Demetrescu 2014, Fan & Lee 2019, AHW 2009, B&D 2015, HLT 2021, LYCP 2019, CY 2006, Elliott 2011) | `ivx_ra(horizon)`, `ivx(lag_y)`, `ivx_qr_boot()`, `arm()`, `ivx_iv()`, `hlt_test()`, `el_test()`, `cy_test()`, `elliott_cf()` | within MC error of the papers; B&D power split differs; see folder README |

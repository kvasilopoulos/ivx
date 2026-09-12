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

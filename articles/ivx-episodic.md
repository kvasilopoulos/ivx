# Episodic predictability: subsample IVX tests

``` r

library(ivx)
```

Full-sample tests find little return predictability; a growing
literature argues that predictability comes and goes in “pockets”.
[`ivx_episodic()`](https://kvasilopoulos.github.io/ivx/reference/ivx_episodic.md)
implements the subsample IVX tests of Demetrescu, Georgiev, Rodrigues &
Taylor (2023, Section 3.2), which formalise the rolling and recursive
approaches of Demetrescu et al. (2022) and, for one-sided tests, of
Pavlidis, Paya & Peel (2017) (see
[`vignette("rolling-ivx")`](https://kvasilopoulos.github.io/ivx/articles/rolling-ivx.md)
for the latter).

## Subsample statistic

For a window \\t = \lfloor \tau_1 T \rfloor + 1, \dots, \lfloor \tau_2 T
\rfloor\\ the IVX statistic is computed from the window’s observations
but with the **full-sample** instrument (DGRT eqs 15–17):

\\ \hat\beta\_{zx}(\tau_1, \tau_2) = \frac{\sum\_{t} z\_{t-1}\\(y_t -
\bar y(\tau_1,\tau_2))} {\sum\_{t} z\_{t-1}\\(x\_{t-1} - \bar
x\_{-1}(\tau_1,\tau_2))}, \qquad t\_{zx}(\tau_1, \tau_2) =
\frac{\hat\beta\_{zx}(\tau_1,\tau_2)}{\mathrm{s.e.}(\hat\beta\_{zx}(\tau_1,\tau_2))},
\\

with \\\hat\sigma_u^2\\ from the window’s OLS residuals (or Eicker–White
weights with `robust = TRUE`). Three agnostic sequences are considered:

- **forward recursive**, \\\\t\_{zx}(0, \tau)\\\\ for \\\tau \in
  \[\tau_L, 1\]\\ — pockets that start at the beginning of the sample;
- **backward recursive**, \\\\t\_{zx}(\tau, 1)\\\\ for \\\tau \in \[0,
  \tau_U\]\\ — end-of-sample pockets;
- **rolling**, \\\\t\_{zx}(\tau, \tau + \Delta\tau)\\\\ for a fixed
  window fraction \\\Delta\tau\\.

The tests are the maximum (right-tailed, \\H_1: \beta \> 0\\), minimum
(left-tailed) and maximum squared (two-sided) of the sequence; with
several predictors, the maximum of the subsample Wald statistics (Remark
11). Their limits are functionals of Brownian motions, so critical
values come from the wild bootstrap of
[`vignette("robust-inference")`](https://kvasilopoulos.github.io/ivx/articles/robust-inference.md):
the same functional is computed on each bootstrap sample. The default is
the fixed regressor wild bootstrap used by Demetrescu et al. (2022);
`type = "rwb"` is available.

``` r

mod <- ivx(Ret ~ DP, data = kms)
e <- ivx_episodic(mod, scheme = "rolling", window = 0.2, B = 499, seed = 1)
e
#> 
#> Call:
#> ivx(formula = Ret ~ DP, data = kms, horizon = 1)
#> 
#> Subsample IVX tests, rolling scheme (window = 0.2), 827 windows
#> Fixed regressor wild bootstrap, B = 499
#> 
#>                         statistic bootstrap p
#> sup t   (H1: beta > 0)      2.921      0.1363
#> inf t   (H1: beta < 0)    -0.4166      0.9980
#> sup t^2 (H1: beta != 0)     8.533      0.2846
```

The sequence itself is returned, with the window bounds as row indices
of the regression sample:

``` r

head(e$sequence, 3)
#>   start end      Wald      t_DP
#> 1     2 207 0.3000535 0.5477714
#> 2     3 208 0.2980866 0.5459731
#> 3     4 209 0.3001673 0.5478752
plot(kms$Date[e$sequence$end], e$sequence$t_DP, type = "l",
     xlab = "window end", ylab = "subsample IVX t-ratio", main = "Rolling windows, 20% of the sample")
abline(h = c(-1.96, 1.96), lty = 2)
```

![](ivx-episodic_files/figure-html/unnamed-chunk-3-1.png)

Forward and backward recursive sequences:

``` r

ivx_episodic(mod, scheme = "forward", window = 0.2, B = 199, seed = 1)$p.value
#>       sup       inf    sup_sq 
#> 0.2462312 0.7788945 0.4371859
ivx_episodic(mod, scheme = "backward", window = 0.8, B = 199, seed = 1)$p.value
#>        sup        inf     sup_sq 
#> 0.03015075 0.98994975 0.05025126
```

## Interpretation and caveats

- A rejection says that predictability of the given sign exists in at
  least one window; the plot of the sequence locates it. The window(s)
  above the pointwise \\\pm 1.96\\ line are not individually significant
  at 5% — that ignores the multiplicity the sup-test accounts for.
- `window` is a fraction of the regression sample; the paper uses window
  fractions between 0.1 and 0.3 for rolling tests and warm-in fractions
  of 0.1–0.25 for recursive ones. Very short windows give noisy
  statistics and low power.
- The instrument is built once from the full sample. Refitting
  [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) on
  each window (as in
  [`vignette("rolling-ivx")`](https://kvasilopoulos.github.io/ivx/articles/rolling-ivx.md))
  rebuilds the instrument inside the window, which is a different
  statistic; the Bonferroni approach used there is conservative, the
  bootstrap sup-test is not.
- Short horizon only; the fit must be a plain `ivx` object without
  weights.

## References

- Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
  (2022). Testing for episodic predictability in stock returns. *Journal
  of Econometrics*, 227(1), 85–113.
- Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
  (2023). Extensions to IVX methods of inference for return
  predictability. *Journal of Econometrics*, 237(2), 105271.
- Pavlidis, E. G., Paya, I., & Peel, D. A. (2017). Testing for
  speculative bubbles using spot and forward prices. *International
  Economic Review*, 58(4), 1191–1226.

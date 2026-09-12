# Rolling IVX tests for bubble detection

Pavlidis, Paya & Peel (2017) detect periodically collapsing bubbles by
running the Fama (1984) regression inside a rolling window and drawing
inference with IVX. Everything the procedure needs is already in
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md), so this
vignette shows how to build it in a few lines rather than through a
dedicated function.

The procedure is:

1.  Fix a window of `w` observations. The paper uses
    `w = max(25, r0 * T)`, where `r0 * T` is the minimum window rule of
    Phillips, Shi & Yu (2015): `r0 = 0.01 + 1.8 / sqrt(T)`.
2.  For every window ending at `t = w, ..., T`, fit the predictive
    regression with
    [`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) and
    keep the IVX t-statistic of the slope. Under the null the statistic
    is standard normal.
3.  Compare the sequence against a constant critical value. For the
    overall test of “any bubble in the sample” the paper applies a
    Bonferroni correction for the `T - w + 1` hypotheses,
    `qnorm(1 - alpha / (T - w + 1))`. For dating the episodes it
    compares each statistic with the plain `qnorm(1 - alpha)`.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ivx`](https://kvasilopoulos.github.io/ivx/)`)`\
\
`ivx_roll`` ``<-`` ``function``(``formula``, ``data``, ``window``, ``horizon`` ``=`` ``1``, ``alpha`` ``=`` ``0.05``, ``...``)`` ``{`\
`  ``n`` ``<-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``data``)`\
`  ``ends`` ``<-`` ``window``:``n`\
`  ``tstat`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``ends``, ``function``(``e``)`` ``{`\
`    ``fit`` ``<-`` `[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``formula``, ``data``[``(``e`` ``-`` ``window`` ``+`` ``1``)``:``e``, ``]``, horizon ``=`` ``horizon``, ``...``)`\
`    ``fit``$``tstat``[``1``]`\
`  ``}``, `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``1``)``)`\
`  `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`    end ``=`` ``ends``,`\
`    tstat ``=`` ``tstat``,`\
`    cv ``=`` `[`qnorm`](https://rdrr.io/r/stats/Normal.html)`(``1`` ``-`` ``alpha``)``,`\
`    cv_bonf ``=`` `[`qnorm`](https://rdrr.io/r/stats/Normal.html)`(``1`` ``-`` ``alpha`` ``/`` `[`length`](https://rdrr.io/r/base/length.html)`(``ends``)``)`\
`  ``)`\
`}`

`fit$tstat` is `coef / se`, i.e. the IVX t-ratio of the first regressor;
the one-sided alternative in the paper is a slope above its
efficient-market value. With the KMS monthly data:

\
`n`` ``<-`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``kms``)`\
`w`` ``<-`` `[`max`](https://rdrr.io/r/base/Extremes.html)`(``25``, `[`floor`](https://rdrr.io/r/base/Round.html)`(``(``0.01`` ``+`` ``1.8`` ``/`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``n``)``)`` ``*`` ``n``)``)`\
`r`` ``<-`` ``ivx_roll``(``Ret`` ``~`` ``LTY``, ``kms``, window ``=`` ``w``)`\
\
[`max`](https://rdrr.io/r/base/Extremes.html)`(``r``$``tstat``)`\
`#> [1] 2.327991`\
`r``$``cv_bonf``[``1``]`\
`#> [1] 3.882191`

The overall test rejects if the maximum statistic exceeds the Bonferroni
critical value. Episodes are dated by the periods where the statistic
sits above the standard normal critical value:

\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``r``$``end``, ``r``$``tstat``, type ``=`` ``"l"``, xlab ``=`` ``"window end"``, ylab ``=`` ``"IVX t-statistic"``)`\
[`abline`](https://rdrr.io/r/graphics/abline.html)`(``h ``=`` ``r``$``cv``[``1``]``, col ``=`` ``"red"``, lty ``=`` ``2``)`\
[`abline`](https://rdrr.io/r/graphics/abline.html)`(``h ``=`` ``r``$``cv_bonf``[``1``]``, col ``=`` ``"red"``)`

![](rolling-ivx_files/figure-html/unnamed-chunk-4-1.png)

\
[`range`](https://rdrr.io/r/base/range.html)`(``r``$``end``[``r``$``tstat`` ``>`` ``r``$``cv``]``)`\
`#> [1] 203 986`

For a long-horizon Fama regression pass `horizon = n` and the
overlapping regressand;
[`ivx()`](https://kvasilopoulos.github.io/ivx/reference/ivx.md) then
uses the Phillips & Lee (2013) long-horizon estimator. Bootstrap
critical values for each window can be obtained by replacing the
[`qnorm()`](https://rdrr.io/r/stats/Normal.html) constants with
[`ivx_boot()`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)
quantiles, at the cost of `B` extra fits per window.

## References

Pavlidis, E. G., Paya, I., & Peel, D. A. (2017). Testing for speculative
bubbles using spot and forward prices. *International Economic Review*,
58(4), 1191-1226.

Phillips, P. C. B., Shi, S., & Yu, J. (2015). Testing for multiple
bubbles: Historical episodes of exuberance and collapse in the S&P 500.
*International Economic Review*, 56(4), 1043-1078.

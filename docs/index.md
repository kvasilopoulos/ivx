# ivx: Robust Econometric Inference

Drawing statistical inference on the coefficients of a short- or
long-horizon predictive regression with persistent regressors by using
the IVX method of [Magdalinos and Phillips
(2009)](doi:10.1017/S0266466608090154) and [Kostakis, Magdalinos and
Stamatogiannis (2015)](doi:10.1093/rfs/hhu139).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

\
`# Install release version from CRAN`\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"ivx"``)`\
\
\
`# install.packages("devtools")`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"kvasilopoulos/ivx"``)`

## Usage

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ivx`](https://kvasilopoulos.github.io/ivx/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`magrittr`](https://magrittr.tidyverse.org)`)`

This is a basic example, lets load the data first:

\
`# Monthly data from Kostakis et al (2014)`\
`kms`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`names`](https://rdrr.io/r/base/names.html)`(``)`\
`#>  [1] "Date" "DE"   "LTY"  "DY"   "DP"   "TBL"  "EP"   "BM"   "INF"  "DFY" `\
`#> [11] "NTIS" "TMS"  "Ret"`

## Univariate

And then do the univariate estimation:

\
[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``Ret`` ``~`` ``DP``, data ``=`` ``kms``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(``)`\
`#> `\
`#> Call:`\
`#> ivx(formula = Ret ~ DP, data = kms, horizon = 1)`\
`#> `\
`#> Coefficients:`\
`#>    Estimate Std. Error t value Wald Ind Pr(> chi)`\
`#> DP 0.006489   0.004553   1.425    2.031     0.154`\
`#> `\
`#> Joint Wald statistic:  2.031 on 1 DF, p-value 0.1541`\
`#> Multiple R-squared:  0.002844,   Adjusted R-squared:  0.001877`\
\
[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``Ret`` ``~`` ``DP``, data ``=`` ``kms``, horizon ``=`` ``4``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(``)`\
`#> `\
`#> Call:`\
`#> ivx(formula = Ret ~ DP, data = kms, horizon = 4)`\
`#> `\
`#> Coefficients:`\
`#>    Estimate Std. Error t value Wald Ind Pr(> chi)`\
`#> DP 0.006931   0.004599   1.507    2.271     0.132`\
`#> `\
`#> Joint Wald statistic:  2.271 on 1 DF, p-value 0.1318`\
`#> Multiple R-squared:  0.01167,    Adjusted R-squared:  0.01358`

## Multivariate

And the multivariate estimation, for one or multiple horizons:

\
[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``Ret`` ``~`` ``DP`` ``+`` ``TBL``, data ``=`` ``kms``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(``)`\
`#> `\
`#> Call:`\
`#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 1)`\
`#> `\
`#> Coefficients:`\
`#>      Estimate Std. Error t value Wald Ind Pr(> chi)`\
`#> DP   0.006145   0.004557   1.349    1.819     0.177`\
`#> TBL -0.080717   0.057701  -1.399    1.957     0.162`\
`#> `\
`#> Joint Wald statistic:  3.644 on 2 DF, p-value 0.1617`\
`#> Multiple R-squared:  0.004968,   Adjusted R-squared:  0.003036`\
\
[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``Ret`` ``~`` ``DP`` ``+`` ``TBL``, data ``=`` ``kms``, horizon ``=`` ``4``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(``)`\
`#> `\
`#> Call:`\
`#> ivx(formula = Ret ~ DP + TBL, data = kms, horizon = 4)`\
`#> `\
`#> Coefficients:`\
`#>      Estimate Std. Error t value Wald Ind Pr(> chi)`\
`#> DP   0.006579   0.004601   1.430    2.045     0.153`\
`#> TBL -0.073549   0.058238  -1.263    1.595     0.207`\
`#> `\
`#> Joint Wald statistic:  3.527 on 2 DF, p-value 0.1715`\
`#> Multiple R-squared:  0.018,  Adjusted R-squared:  0.01895`

## Robust inference (Demetrescu et al., 2023)

Eicker-White standard errors, the IVX tuning parameters and
wild-bootstrap p-values (residual or fixed-regressor wild bootstrap):

\
[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``Ret`` ``~`` ``DP`` ``+`` ``TBL``, data ``=`` ``kms``, robust ``=`` ``TRUE``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(``)`\
`#> `\
`#> Call:`\
`#> ivx(formula = Ret ~ DP + TBL, data = kms, robust = TRUE, horizon = 1)`\
`#> `\
`#> Coefficients:`\
`#>      Estimate Std. Error t value Wald Ind Pr(> chi)`\
`#> DP   0.006145   0.004792   1.282    1.644     0.200`\
`#> TBL -0.080717   0.057250  -1.410    1.988     0.159`\
`#> (Eicker-White standard errors)`\
`#> `\
`#> Joint Wald statistic:  2.893 on 2 DF, p-value 0.2354`\
`#> Multiple R-squared:  0.004968,   Adjusted R-squared:  0.003036`\
\
`mod`` ``<-`` `[`ivx`](https://kvasilopoulos.github.io/ivx/reference/ivx.md)`(``Ret`` ``~`` ``DP`` ``+`` ``TBL``, data ``=`` ``kms``, beta ``=`` ``0.9``, cz ``=`` ``5``)`\
[`ivx_boot`](https://kvasilopoulos.github.io/ivx/reference/ivx_boot.md)`(``mod``, B ``=`` ``999``, type ``=`` ``"rwb"``, seed ``=`` ``1``)`\
`#> `\
`#> Call:`\
`#> ivx(formula = Ret ~ DP + TBL, data = kms, beta = 0.9, cz = 5, `\
`#>     horizon = 1)`\
`#> `\
`#> Residual wild bootstrap, B = 999`\
`#> `\
`#> Coefficients (bootstrap p-values):`\
`#>     Estimate t value Wald Ind Pr(> chi) Pr(t < 0) Pr(t > 0)`\
`#> DP   0.00502   0.925    0.855    0.5415   0.60761     0.392`\
`#> TBL -0.12896  -1.581    2.498    0.1652   0.08308     0.917`\
`#> `\
`#> Joint Wald statistic: 2.807, bootstrap p-value 0.4114`

## Yang et al. (2020) IVX-AR methodology

\
[`ivx_ar`](https://kvasilopoulos.github.io/ivx/reference/ivx_ar.md)`(``hpi`` ``~`` ``cpi``, data ``=`` ``ylpc``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `\
`  `[`summary`](https://rdrr.io/r/base/summary.html)`(``)`\
`#> `\
`#> Call:`\
`#> ivx_ar(formula = hpi ~ cpi, data = ylpc, horizon = 1)`\
`#> `\
`#> Auto () with AR terms q = 4`\
`#> `\
`#> Coefficients:`\
`#>       Estimate Std. Error t value Wald Ind Pr(> chi)  `\
`#> cpi -1.775e-04  8.532e-05  -2.080    4.326    0.0375 *`\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> Joint Wald statistic:  4.326 on 1 DF, p-value 0.03753`\
`#> Multiple R-squared:  0.02721,    Adjusted R-squared:  0.02142`\
`#> Wald AR statistic: 132.3 on 4 DF, p-value < 2.2e-16`

------------------------------------------------------------------------

Please note that the ‘ivx’ project is released with a [Contributor Code
of
Conduct](https://github.com/kvasilopoulos/ivx/blob/master/.github/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

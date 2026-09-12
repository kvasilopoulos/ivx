# Moving Block Bootstrap for IVX-QR

Percentile confidence intervals and p-values for the IVX-QR coefficients
from the moving block bootstrap (MBB) of Fan and Lee (2019, Section 5).
Blocks of the pairs \\(y_t, ilde z\_{t-1})\\ are resampled and the
quantile regression of Lee (2016) is refitted on each bootstrap sample.
This avoids estimating the sparsity and the nuisance parameters that
appear under conditional heteroskedasticity, which is where the
asymptotic IVX-QR test is most distorted, especially in the tails.

## Usage

``` r
ivx_qr_boot(object, B = 999, block = NULL, level = 0.95, seed = NULL)

# S3 method for class 'ivx_qr_boot'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  an object of class "ivx_qr".

- B:

  number of bootstrap replications.

- block:

  block length; the default is \\\lceil n^{1/4} \rceil\\ as in the
  paper.

- level:

  confidence level of the percentile intervals.

- seed:

  optional integer seed.

- x:

  an object of class "ivx_qr_boot".

- digits:

  minimal number of significant digits.

- ...:

  unused.

## Value

an object of class "ivx_qr_boot": a list with the estimates, the
percentile intervals `ci`, two-sided percentile p-values `p.value` for
\\H_0: eta_j = 0\\, and the bootstrap draws `boot`.

## References

Fan, R., & Lee, J. H. (2019). Predictive quantile regressions under
persistence and conditional heteroskedasticity. Journal of Econometrics,
213(1), 261-280.

## Examples

``` r
if (requireNamespace("quantreg", quietly = TRUE)) {
  m <- ivx_qr(Ret ~ DP, data = kms, tau = 0.1)
  ivx_qr_boot(m, B = 199, seed = 1)
}
#> 
#> Call:
#> ivx_qr(formula = Ret ~ DP, data = kms, tau = 0.1)
#> 
#> IVX-QR at tau = 0.1, moving block bootstrap, B = 199, block length 6
#> 
#> Coefficients (percentile intervals and p-values):
#>     Estimate      2.5%     97.5% Pr(|b| > 0)
#> DP -0.024136 -0.057003  0.006266       0.131
#> 
```

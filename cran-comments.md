## Release summary

Update of ivx 1.1.1 (on CRAN) to 2.0.0. New methods (residual-augmented,
lag-augmented and long-horizon IVX, IVX quantile regression with a block
bootstrap, systems IVX, wild-bootstrap and subsample inference, combined
instruments, and five non-IVX benchmark tests), bug fixes in `ac_test_bg()`,
`ac_test_lb()`/`ac_test_bp()` and weighted fits, and eight vignettes. See
NEWS.md.

`drop1()`, `add1()` and `step()` methods for `ivx` objects are deprecated
(they still work, with a warning). No function was removed or renamed.

## Test environments

* local Windows 11, R 4.6.1 (R CMD check --as-cran on the built tarball)
* GitHub Actions: ubuntu-latest (R release, devel), windows-latest (release),
  macOS-latest (release)
* win-builder: see below

## R CMD check results

0 errors | 0 warnings | 0 notes

(The local Windows run shows one NOTE, "non-standard things in the check
directory: 'NULL'", an artefact of the local check environment; it does not
appear on the CI platforms.)

## Downstream dependencies

There are no reverse dependencies on CRAN.

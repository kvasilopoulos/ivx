# CRAN release checklist — ivx 2.0.0 (update of 1.1.1)

Compiled from the CRAN Repository Policy
(<https://cran.r-project.org/web/packages/policies.html>), the CRAN
submission checklist
(<https://cran.r-project.org/web/packages/submission_checklist.html>),
*R Packages* ch. “Releasing to CRAN” (<https://r-pkgs.org/release.html>)
and `usethis::use_release_issue()`. Tick boxes record the state of the
checks run locally on 2026-09-13; `[ ]` items are for the maintainer.

## Policy points that apply to this package

- Update of an existing CRAN package (1.1.1 → 2.0.0, major). CRAN asks
  for updates “no more than every 1–2 months”; last release was long
  ago.
- 0 reverse dependencies on CRAN
  (`tools::package_dependencies("ivx", reverse = TRUE)`), so no
  `revdepcheck` and no two-week notice to downstream maintainers.
- License GPL-3 unchanged. Maintainer unchanged (Kostas Vasilopoulos).
- Compiled code (Rcpp/RcppArmadillo): no `abort`/`exit`/`assert`; errors
  go through `Rcpp::stop()` only.
- Multi-core: `ivx_boot(cores = )` and `ivx_episodic(cores = )` default
  to 1; tests and examples use ≤ 2 cores (the parallel test uses 2 and
  is `skip_on_cran()`).
- Nothing is written outside
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html); no internet
  access at run time.
- Tarball 0.6 MB; data 5 files, \< 1 MB; vignettes HTML.
- Examples: all under 5 s each (`--timings`); the bootstrap examples use
  `B = 199`.
- API change:
  [`drop1()`](https://rdrr.io/r/stats/add1.html)/[`add1()`](https://rdrr.io/r/stats/add1.html)/[`step()`](https://rdrr.io/r/stats/step.html)
  on `ivx` objects are deprecated (warn, still work). No function was
  removed or renamed.

## Pre-submission checks

`devtools::document()`; NAMESPACE and Rd regenerated

`devtools::test()` — 0 failures (parallel bootstrap test skipped on
CRAN)

[`covr::package_coverage()`](http://covr.r-lib.org/reference/package_coverage.md)
— 94%

[`spelling::spell_check_package()`](https://docs.ropensci.org/spelling//reference/spell_check_package.html)
— 0 (Language en-GB, `inst/WORDLIST` updated)

`urlchecker::url_check()` — all URLs OK. The DOI badge in README now
links to `doi.org` (concept DOI 10.5281/zenodo.3371391) rather than
`zenodo.org`, which timed out in earlier checks (HTTP 504 in the CRAN
incoming check too); only the badge image still comes from zenodo.org.

After tagging the GitHub release, confirm Zenodo archived it (the
GitHub-Zenodo hook did not fire for v1.1.1): check
<https://zenodo.org/account/settings/github/> and the concept record.
`.zenodo.json` supplies the deposit metadata; `CITATION.cff` needs
`version` and `date-released` bumped with each release.

`R CMD build` with vignettes (pandoc from RStudio, TinyTeX for the
manual)

`R CMD check --as-cran` on the tarball with current R release: see
`cran-comments.md` for the result of the last run

`devtools::check_win_devel()` and `check_win_release()` — uploads to
win-builder and emails the maintainer (not run automatically; run it
yourself, or say so and it will be launched)

R-hub v2 (`rhub::rhub_setup()` + `rhub::rhub_check()` on GitHub Actions)
for Linux/macOS/R-devel — optional; the GitHub Actions matrix already
covers ubuntu (release, devel), windows, macOS

Check the GitHub Actions run on the release commit is green

`NEWS.md`: 2.0.0 section complete (it is; re-read once more)

`DESCRIPTION`: Version 2.0.0, `Date` field optional (omitted)

`cran-comments.md` updated (done below)

`README.md` re-knit from `README.Rmd` (done; re-knit if anything
changes)

## Submission

1.  `devtools::submit_cran()` (builds the tarball, uploads it with
    `cran-comments.md`), or the web form
    <https://cran.r-project.org/submit.html>.
2.  Confirm the e-mail CRAN sends to the maintainer address.
3.  Watch <https://cran.r-project.org/incoming/> and the maintainer
    inbox.
4.  If CRAN asks for changes: fix, bump to 1.2.1, add a “Resubmission”
    section at the top of `cran-comments.md` describing what changed,
    resubmit.

## After acceptance

- `usethis::use_github_release()` (release notes from NEWS.md), tag
  `v2.0.0`
- `usethis::use_dev_version(push = TRUE)` → 2.0.0.9000
- Wait 48 h for the CRAN check page before submitting any correction

## Files that must NOT ship

Excluded via `.Rbuildignore`: `research/`, `data-raw/`, `docs/`,
`pkgdown/`, `_pkgdown.yml`, `README.Rmd`, `CLAUDE.md`,
`cran-release.md`, `cran-comments.md`, `CITATION.cff`, `.github/`,
`.env`, gcov files, `Rplots.pdf`, `codecov.yml`. Verify with
`tar tzf ivx_2.0.0.tar.gz`.

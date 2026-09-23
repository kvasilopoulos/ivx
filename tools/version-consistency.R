#!/usr/bin/env Rscript
# Check that the package version is consistent across the files that repeat it.
#
# DESCRIPTION is the source of truth. NEWS.md must always agree with it.
# CITATION.cff records the last *released* version, so it is only checked when
# DESCRIPTION holds a release version (dev versions carry a fourth component).
#
# Usage:
#   Rscript tools/version-consistency.R          # check, non-zero exit on drift
#   Rscript tools/version-consistency.R --fix    # also sync CITATION.cff

args <- commandArgs(trailingOnly = TRUE)
fix <- "--fix" %in% args

desc_version <- unname(read.dcf("DESCRIPTION", fields = "Version")[1, 1])
is_dev <- length(unlist(strsplit(desc_version, ".", fixed = TRUE))) > 3

news <- readLines("NEWS.md", warn = FALSE)
news_heading <- grep("^#[[:space:]]+ivx[[:space:]]+", news, value = TRUE)[1]
news_version <- sub("^#[[:space:]]+ivx[[:space:]]+", "", news_heading)

cff <- readLines("CITATION.cff", warn = FALSE)
cff_line <- grep("^version:", cff)[1]
cff_version <- trimws(sub("^version:", "", cff[cff_line]))

problems <- character()

if (is.na(news_version) || !identical(news_version, desc_version)) {
  problems <- c(problems, sprintf(
    "NEWS.md top heading is '%s' but DESCRIPTION Version is '%s'",
    news_heading, desc_version
  ))
}

if (is_dev) {
  message(sprintf(
    "DESCRIPTION is a development version (%s); CITATION.cff (%s) not checked.",
    desc_version, cff_version
  ))
} else if (!identical(cff_version, desc_version)) {
  if (fix) {
    cff[cff_line] <- paste0("version: ", desc_version)
    date_line <- grep("^date-released:", cff)[1]
    cff[date_line] <- paste0("date-released: ", format(Sys.Date()))
    writeLines(cff, "CITATION.cff")
    message(sprintf(
      "CITATION.cff: version %s -> %s, date-released -> %s",
      cff_version, desc_version, format(Sys.Date())
    ))
  } else {
    problems <- c(problems, sprintf(
      "CITATION.cff version is '%s' but DESCRIPTION Version is '%s' (run with --fix)",
      cff_version, desc_version
    ))
  }
}

if (length(problems)) {
  message("Version drift:\n", paste0("  - ", problems, collapse = "\n"))
  quit(status = 1L)
}

message(sprintf(
  "Version %s consistent across DESCRIPTION, NEWS.md and CITATION.cff.",
  desc_version
))

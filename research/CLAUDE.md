# research/

Literature for package extensions. Not part of the R build (`.Rbuildignore`).

- `README.md` — the index. Paper table (file basename, version, reference, why it matters),
  "Candidates not yet collected" (ranked, with DOIs), "Suggested extensions".
- `pdf/<basename>.pdf` and `txt/<basename>.txt` share basenames; `txt/` is
  `pdftotext -layout pdf/X.pdf txt/X.txt`. Grep `txt/` when checking a formula or algorithm.
- Naming: `author(s)-year-short-slug`; arXiv-only papers keep `arxiv-<id>-slug`.
- Prefer the published journal version; the "version" column says what is in `pdf/`.

Adding a paper: drop the PDF in `pdf/`, rename to the convention, run `pdftotext -layout`,
add a table row (or move it out of "Candidates not yet collected").

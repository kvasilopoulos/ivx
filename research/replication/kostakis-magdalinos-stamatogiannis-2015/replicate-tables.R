# Replication of Kostakis, Magdalinos & Stamatogiannis (2015, RFS 28, 1506-1553)
# with the `kms` dataset shipped in the package (monthly, 1927-2012).
# The same numbers are asserted in tests/testthat/test-ivx.R.
#
# Usage: Rscript replicate-tables.R

here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
pkg <- normalizePath(file.path(here, "..", "..", ".."))
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(pkg, quiet = TRUE) else library(ivx)

wald <- function(m) unname(coef(summary(m))[, "Wald Ind"])
show <- function(title, paper, ivx) {
  cat("\n", title, "\n", sep = "")
  print(data.frame(paper = paper, ivx = round(ivx, 3), diff = round(ivx, 3) - paper), row.names = FALSE)
}

# Table 6 (p. 1531): univariate IVX, Ret ~ DE
m <- ivx(Ret ~ DE, data = kms)
show("Table 6, Ret ~ DE: coefficient, Wald, delta",
     c(-0.0033, 0.393, -0.067), c(round(coef(m), 4), wald(m), delta(m)))

# Table 8 (p. 1537): multivariate IVX, Ret ~ DP + TBL
m <- ivx(Ret ~ DP + TBL, data = kms)
show("Table 8, Ret ~ DP + TBL: coefficients, joint Wald",
     c(0.0061, -0.0807, 3.644), c(round(coef(m), 4), m$Wald_Joint))

# Table 11 (p. 1544): long-horizon univariate Wald, Ret ~ DE
h <- c(4, 12, 24, 36, 48, 60)
show("Table 11, Ret ~ DE, Wald by horizon (4, 12, 24, 36, 48, 60)",
     c(0.138, 0.005, 0.472, 0.803, 0.422, 0.637),
     sapply(h, function(k) wald(ivx(Ret ~ DE, data = kms, horizon = k))))

# Table 13 (p. 1547): long-horizon multivariate, Ret ~ EP + TBL
tab13 <- t(sapply(h, function(k) { m <- ivx(Ret ~ EP + TBL, data = kms, horizon = k); c(wald(m), m$Wald_Joint) }))
show("Table 13, Ret ~ EP + TBL: Wald EP, Wald TBL, joint Wald by horizon (rows = 4, 12, 24, 36, 48, 60)",
     c(5.778, 6.383, 4.990, 4.599, 4.983, 4.321, 3.894, 3.166, 2.124, 1.915, 1.441, 1.039,
       7.638, 7.614, 5.794, 5.383, 5.660, 4.822),
     as.vector(tab13))

# Monte Carlo size check of ivx_ra() against Demetrescu & Rodrigues (2022, JoE 227,
# 429-460), Table 3 (right-sided tests, b = 0 rows) and Table 4 (two-sided), T = 200.
# DGP (24)-(25): y_t = beta x_{t-1} + u_t, x_t = rho x_{t-1} + v_t, v_t = a1 v_{t-1} + e_t,
# rho = 1 - c/T, a1 = -0.5, corr(u_t, e_t) = -0.95, nominal level 5%.
#
# Usage: Rscript size-table3.R [R]   (paper: 10000 replications; default 2000, ~2 min)

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args)) as.integer(args[1]) else 2000L
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
pkg <- normalizePath(file.path(here, "..", "..", ".."))
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(pkg, quiet = TRUE) else library(ivx)

sim <- function(Tn, c, a1 = -0.5, r = -0.95) {
  e <- rnorm(Tn)
  u <- r * e + sqrt(1 - r^2) * rnorm(Tn)
  v <- as.numeric(stats::filter(e, a1, "recursive"))
  x <- as.numeric(stats::filter(v, 1 - c / Tn, "recursive"))
  list(y = c(0, u[-1]), x = matrix(x))
}
paper <- data.frame(
  c = c(0, 10, 20, 30, 40, 50),
  right_ivx = c(0.116, 0.088, 0.074, 0.066, 0.064, 0.061),
  right_ra  = c(0.054, 0.055, 0.055, 0.053, 0.050, 0.050)
)
set.seed(2022)
Tn <- 200
res <- t(sapply(paper$c, function(cc) {
  tt <- replicate(R, {
    d <- sim(Tn, cc)
    c(ivx_fit(d$y, d$x)$tstat, ivx_ra_fit(d$y, d$x)$tstat)
  })
  c(right_ivx = mean(tt[1, ] > qnorm(0.95)), right_ra = mean(tt[2, ] > qnorm(0.95)),
    two_ivx = mean(abs(tt[1, ]) > qnorm(0.975)), two_ra = mean(abs(tt[2, ]) > qnorm(0.975)))
}))
out <- cbind(paper, ivx = round(res, 3))
print(out, row.names = FALSE)
cat("\nMC standard error at 5% with R =", R, ":", round(sqrt(0.05 * 0.95 / R), 3), "\n")
write.csv(out, file.path(here, "size-table3.csv"), row.names = FALSE)

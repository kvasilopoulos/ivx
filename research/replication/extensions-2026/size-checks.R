# Monte Carlo checks for the methods added in ivx 1.2.0, reproducing the numbers quoted
# in the vignettes (ivx.Rmd, ivx-ra.Rmd, ivx-qr.Rmd, ivx-iv.Rmd, benchmarks.Rmd). Each block
# names the paper table or statement it is compared with. Nominal level 5% unless noted.
#
# Usage: Rscript size-checks.R [R]   (default 500 replications per cell, ~15 min)

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args)) as.integer(args[1]) else 500L
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
pkg <- normalizePath(file.path(here, "..", "..", ".."))
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(pkg, quiet = TRUE) else library(ivx)
out <- list()
rej <- function(expr) { e <- substitute(expr); env <- parent.frame(); mean(replicate(R, eval(e, env))) }
dgp <- function(n, c = 0, d = -0.95, b = 0, x0 = 0) {
  e <- rnorm(n); u <- d * e + sqrt(1 - d^2) * rnorm(n)
  x <- as.numeric(stats::filter(e, 1 + c / n, "recursive", init = x0))
  list(y = c(0, b / n * x[-n] + u[-1]), x = x, e = e, u = u)
}

# 1. ivx_ra(horizon): DRT (2023) report two-sided sizes in [0.023, 0.058] for T = 500, delta = -0.95
set.seed(1)
out$ivx_ra_horizon <- sapply(c(0, -20), function(cc) sapply(c(10, 20), function(h)
  rej({ s <- dgp(500, cc); abs(ivx_ra_fit(s$y, matrix(s$x), ar = 1, horizon = h)$tstat) > 1.96 })))
dimnames(out$ivx_ra_horizon) <- list(h = c(10, 20), c = c(0, -20))

# 2. ivx(lag_y): Demetrescu (2014) Table 1, T = 100, delta = -0.98, eta = 0 (paper: 4.6/7.9 at b = 0, 8.3/63.9 at b = 5)
set.seed(3)
out$lag_y <- sapply(c(0, 5), function(b) c(
  ivx = rej({ s <- dgp(100, 0, -0.98, b); abs(ivx_fit(s$y, matrix(s$x), beta = 0)$tstat[1]) > 1.96 }),
  lag_y = rej({ s <- dgp(100, 0, -0.98, b); abs(ivx_fit(s$y, matrix(s$x), beta = 0, lag_y = TRUE)$tstat[1]) > 1.96 })))
colnames(out$lag_y) <- c("b0", "b5")

# 3. ivx_qr_boot: Fan & Lee (2019) study 2, ARCH(1) errors, tau = 0.1, n = 200 (asymptotic test oversized, MBB near nominal)
if (requireNamespace("quantreg", quietly = TRUE)) {
  set.seed(5)
  out$qr_boot <- rowMeans(replicate(R, {
    n <- 200; e <- rnorm(n); ux <- -0.9 * e + sqrt(0.19) * rnorm(n)
    u <- numeric(n); for (t in 2:n) u[t] <- sqrt(1 + 0.9 * u[t - 1]^2) * e[t]
    x <- cumsum(ux); y <- u - quantile(u, 0.1)
    f <- ivx_qr_fit(y, matrix(x), tau = 0.1); f$cnames <- "x"; class(f) <- c("ivx_qr", "ivx")
    c(asymptotic = f$Wald_Joint > qchisq(0.95, 1), mbb = ivx_qr_boot(f, B = 199)$p.value < 0.05)
  }))
}

# 4. arm: AHW (2009) Case 1, Phi = [.8 .1; .1 .85], Sigma_v = [2 1; 1 2], n = 100 (paper sizes 6-10% at n = 50, 5-7% at n = 200)
set.seed(9)
Phi <- matrix(c(.8, .1, .1, .85), 2); L <- chol(matrix(c(2, 1, 1, 2), 2))
out$arm <- rowMeans(replicate(R, {
  n <- 100; v <- matrix(rnorm(2 * n), n) %*% L
  x <- matrix(0, n, 2); for (t in 2:n) x[t, ] <- Phi %*% x[t - 1, ] + v[t, ]
  y <- -0.9 * rowSums(v) / sqrt(6) + sqrt(0.19) * rnorm(n)
  f <- arm_fit(y, x); c(t1 = abs(f$tstat[1]) > 1.96, t2 = abs(f$tstat[2]) > 1.96, wald = f$Wald_Joint > qchisq(.95, 2))
}))

# 5. ivx_iv: Breitung & Demetrescu (2015) Table 1, T = 250, corr 0.9, 10% level (paper rho=1: comb 11.2/65.7/91.3 at b = 0/10/20)
set.seed(2)
out$ivx_iv <- sapply(c(0, 10, 20), function(b) sapply(c("comb", "sin", "frac", "diff"), function(i)
  rej({ s <- dgp(250, 0, 0.9, b); abs(ivx_iv_fit(s$y, matrix(s$x), i)$tstat) > 1.645 })))
colnames(out$ivx_iv) <- c("b0", "b10", "b20")

# 6. hlt_test: HLT (2021) design, T = 200, upper-tail 5% (paper: near nominal, undersized for positive rho_xy)
set.seed(1)
out$hlt <- sapply(c(-0.9, 0.9), function(r) sapply(c(1, 0.95, 0.5), function(rho)
  rej({ n <- 200; e <- matrix(rnorm(2 * n), n); v <- e[, 1]; u <- r * v + sqrt(1 - r^2) * e[, 2]
        hlt_test_fit(u, as.numeric(stats::filter(v, rho, "recursive", init = rnorm(1))))$reject })))
dimnames(out$hlt) <- list(rho = c(1, 0.95, 0.5), rho_xy = c(-0.9, 0.9))

# 7. el_test: LYCP (2019) Gaussian version, beta1 = 0.5, beta2 = 0, n = 500 (paper: accurate size for beta2)
set.seed(2)
out$el <- sapply(c(1, 0.9), function(rho) rej({
  n <- 500; e <- matrix(rnorm(2 * n), n); v <- e[, 1]; u <- -0.5 * v + sqrt(0.75) * e[, 2]
  x <- as.numeric(stats::filter(v, rho, "recursive")); y <- u; y[3:n] <- 0.5 * (x[2:(n - 1)] - x[1:(n - 2)]) + u[3:n]
  el_test_fit(y, x)$p.value[["beta2"]] < 0.05 }))
names(out$el) <- c("rho1", "rho0.9")

# 8. cy_test: CY (2006), T = 200 (paper Section 3.4: right-tailed ~4-5%, left-tailed as low as 1.2%)
set.seed(5)
out$cy <- sapply(c(-0.9, -0.5), function(d) sapply(c(0, -20), function(cc)
  rowMeans(replicate(R, { s <- dgp(200, cc, d); cy_test_fit(s$y, s$x)$reject }))))
dimnames(out$cy) <- list(c("greater_c0", "less_c0", "greater_c-20", "less_c-20"), delta = c(-0.9, -0.5))

# 9. elliott_cf: orthogonalising covariate z = v + noise, delta = -0.9, c = 0, T = 200
set.seed(8)
out$elliott_cf <- rowMeans(replicate(R, {
  n <- 200; v <- rnorm(n); u <- -0.9 * v + sqrt(0.19) * rnorm(n); z <- v + 0.3 * rnorm(n); x <- cumsum(v)
  c(cf = elliott_cf_fit(u, x, z)$p.value < 0.05, ols = summary(lm(u[-1] ~ x[-n]))$coefficients[2, 4] < 0.05)
}))

sink(file.path(here, "size-checks-output.txt"))
cat("Replications per cell:", R, "\n\n")
for (nm in names(out)) { cat("==", nm, "\n"); print(round(out[[nm]], 3)); cat("\n") }
sink()
saveRDS(out, file.path(here, "size-checks.rds"))
cat("done\n")

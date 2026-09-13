# Local-to-unity null quantiles of the DF-GLS statistic (Elliott, Rothenberg &
# Stock, 1996; c_bar = -7, intercept case) on a grid of c, used by cy_test() to
# invert the statistic into a confidence interval for c as in Stock (1991) and
# Campbell & Yogo (2006). Discretised OU process with n = 600 steps, 20,000
# replications per c. Result: `dfgls_q`, a matrix of quantiles (rows: probabilities
# `p`, columns: `c`), stored in R/sysdata.rda.
set.seed(20260913)
n <- 600
R <- 20000
cs <- c(seq(-100, -40, by = 4), seq(-38, 10, by = 1))
p <- seq(0.005, 0.995, by = 0.005)
rb <- 1 - 7 / n

dfgls_sim <- function(c) unlist(lapply(seq_len(R / 5000), function(i) dfgls_chunk(c, 5000)))

dfgls_chunk <- function(c, R) {
  rho <- 1 + c / n
  e <- matrix(rnorm(n * R), n, R)
  x <- e
  for (t in 2:n) x[t, ] <- rho * x[t - 1, ] + e[t, ]
  # quasi-GLS demeaning: regress (x_1, x_t - rb x_{t-1}) on (1, 1 - rb)
  z1 <- x[1, ]
  zt <- x[-1, ] - rb * x[-n, ]
  mu <- (z1 + (1 - rb) * colSums(zt)) / (1 + (n - 1) * (1 - rb)^2)
  xd <- sweep(x, 2, mu)
  xl <- xd[-n, ]
  dx <- xd[-1, ] - xl
  sxx <- colSums(xl^2)
  theta <- colSums(xl * dx) / sxx
  res <- dx - sweep(xl, 2, theta, `*`)
  s2 <- colSums(res^2) / (n - 2)
  theta / sqrt(s2 / sxx)
}

dfgls_q <- sapply(cs, function(c) {
  cat("c =", c, "\n")
  stats::quantile(dfgls_sim(c), p, names = FALSE)
})
dimnames(dfgls_q) <- list(p = p, c = cs)

# Campbell & Yogo (2006), Table 2: DF-GLS confidence-interval levels (lower a1,
# upper a1) for the 5% one-sided Bonferroni Q-test (alpha_2 = 0.10), by delta
cy_table2 <- cbind(
  delta = -c(0.999, seq(0.975, 0.025, by = -0.025)),
  a1_lower = c(0.050, 0.055, 0.055, 0.055, 0.060, 0.060, 0.060, 0.060, 0.065, 0.065,
               0.065, 0.065, 0.070, 0.070, 0.070, 0.075, 0.075, 0.075, 0.080, 0.080,
               0.080, 0.085, 0.085, 0.090, 0.090, 0.095, 0.100, 0.100, 0.105, 0.110,
               0.115, 0.125, 0.130, 0.140, 0.150, 0.160, 0.175, 0.190, 0.215, 0.250),
  a1_upper = c(0.055, 0.080, 0.100, 0.115, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180,
               0.190, 0.195, 0.205, 0.215, 0.225, 0.230, 0.240, 0.250, 0.260, 0.270,
               0.280, 0.285, 0.295, 0.310, 0.320, 0.330, 0.345, 0.355, 0.360, 0.370,
               0.375, 0.380, 0.390, 0.395, 0.400, 0.405, 0.415, 0.420, 0.425, 0.435)
)
usethis::use_data(dfgls_q, cy_table2, internal = TRUE, overwrite = TRUE)

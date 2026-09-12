# Replication of Demetrescu, Georgiev, Rodrigues & Taylor (2023, JoE 237, 105271),
# Table 4, Panel A: IVX predictability tests for the equity premium with the
# Welch-Goyal (2008) monthly predictors, January 1927 - December 2020 (T = 1128).
#
# Paper conventions (Section 6.1):
#   * asymptotic IVX tests use Eicker-White standard errors  -> ivx(robust = TRUE)
#   * bootstrap IVX tests use the residual wild bootstrap with conventional
#     standard errors                                         -> ivx_boot(type = "rwb")
#   * 9999 bootstrap replications (set B below; 1999 takes ~10 min per return series)
#
# The paper describes the dependent variable as the log S&P 500 return including
# dividends minus the log risk-free rate (CRSP_SPvw in the Goyal-Welch file).
# See README.md: the published p-values for dp, ep and bm are only reproduced
# with the ex-dividend series (CRSP_SPvwx), so both are run and reported.
#
# Usage: Rscript replicate-table4.R [B] [vw|vwx|both]
# Run from anywhere; paths are resolved relative to this file.

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 1999L
ret_series <- if (length(args) >= 2) args[2] else "both"
ret_series <- if (ret_series == "both") c("vw", "vwx") else ret_series

here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
pkg <- normalizePath(file.path(here, "..", "..", ".."))
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(pkg, quiet = TRUE) else library(ivx)

# ---- data: Goyal-Welch monthly sheet (see README for source and vintage) ----
gw <- read.csv(file.path(here, "goyal-welch-monthly.csv"),
               na.strings = c("NaN", ""), check.names = FALSE)
gw$Index <- as.numeric(gsub(",", "", gw$Index))
lag1 <- function(x) c(NA, head(x, -1))
d <- data.frame(
  Date = as.Date(paste0(gw$yyyymm, "01"), "%Y%m%d"),
  vw   = log(1 + gw$CRSP_SPvw) - log(1 + gw$Rfree),   # log excess return, with dividends
  vwx  = log(1 + gw$CRSP_SPvwx) - log(1 + gw$Rfree),  # log excess return, ex dividends
  dp   = log(gw$D12) - log(gw$Index),
  dy   = log(gw$D12) - log(lag1(gw$Index)),
  ep   = log(gw$E12) - log(gw$Index),
  de   = log(gw$D12) - log(gw$E12),
  svar = gw$svar, bm = gw$`b/m`, ntis = gw$ntis, tbl = gw$tbl, lty = gw$lty,
  ltr  = gw$ltr, tms = gw$lty - gw$tbl, dfy = gw$BAA - gw$AAA,
  dfr  = gw$corpr - gw$ltr, infl = gw$infl
)
# ivx() lags the predictor internally, so keep 1926-12 (predictor) .. 2020-12 (return)
d <- d[d$Date >= as.Date("1926-12-01") & d$Date <= as.Date("2020-12-01"), ]
stopifnot(nrow(d) == 1129)

# ---- published values (Table 4, Panel A) ----
paper <- read.csv(text = "
pred,t2sls_rwb,t2sls_ew,rwb_left,rwb_right,rwb_two,ew_left,ew_right,ew_two,b_ols,b_ivx,rho,delta
dp,0.5125,0.4916,0.5366,0.4634,0.6124,0.7450,0.2550,0.5099,0.002,0.003,0.994,-0.980
dy,0.3395,0.2578,0.8355,0.1645,0.2868,0.9254,0.0746,0.1492,0.000,0.000,1.006,-0.060
ep,0.9116,0.9107,0.8271,0.1729,0.2719,0.9142,0.0858,0.1716,0.005,0.006,0.989,-0.767
de,0.2721,0.2684,0.2236,0.7764,0.4647,0.3013,0.6987,0.6025,-0.004,-0.004,0.991,-0.073
svar,0.5812,0.5821,0.3511,0.6489,0.7107,0.3856,0.6144,0.7712,-0.159,-0.189,0.577,-0.297
bm,0.5819,0.5805,0.6720,0.3280,0.4530,0.7181,0.2819,0.5939,0.008,0.007,0.987,-0.821
ntis,0.5083,0.5061,0.0784,0.9216,0.1609,0.1170,0.8830,0.2339,-0.142,-0.141,0.981,-0.047
tbl,0.3212,0.3135,0.0340,0.9660,0.0874,0.0365,0.9635,0.0730,-0.001,-0.001,0.994,-0.053
lty,0.1145,0.1073,0.0278,0.9722,0.0703,0.0272,0.9728,0.0543,-0.001,-0.001,0.997,-0.088
ltr,0.4045,0.3992,0.9176,0.0824,0.1707,0.9180,0.0820,0.1639,0.001,0.001,0.043,0.055
tms,0.8831,0.8855,0.7883,0.2117,0.4144,0.7434,0.2566,0.5132,0.002,0.001,0.962,-0.002
dfy,0.7160,0.7121,0.4588,0.5412,0.9939,0.5019,0.4981,0.9962,0.000,0.000,0.975,-0.265
dfr,0.2367,0.2316,0.8166,0.1834,0.3740,0.8029,0.1971,0.3943,0.002,0.002,-0.102,0.185
infl,0.5870,0.5887,0.0959,0.9041,0.1874,0.0567,0.9433,0.1113,-0.004,-0.005,0.480,0.033
", strip.white = TRUE)

# ---- replication ----
run_one <- function(ret) {
  set.seed(20230101)
  rows <- lapply(paper$pred, function(p) {
    f <- as.formula(paste(ret, "~", p))
    m_ew <- ivx(f, data = d, robust = TRUE)
    m    <- ivx(f, data = d)
    t_ew <- m_ew$tstat
    b <- ivx_boot(m, B = B, type = "rwb")
    data.frame(
      ret = ret, pred = p,
      rwb_left = b$p.value$tstat[, "less"], rwb_right = b$p.value$tstat[, "greater"],
      rwb_two = b$p.value$tstat[, "two-sided"],
      ew_left = pnorm(t_ew), ew_right = 1 - pnorm(t_ew), ew_two = 2 * (1 - pnorm(abs(t_ew))),
      b_ols = m$ols$coefficients[2], b_ivx = coef(m), rho = m$AR$Rn, delta = delta(m),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}
rep <- do.call(rbind, lapply(ret_series, run_one))
rep$B <- B
write.csv(rep, file.path(here, "table4-panelA-ivx.csv"), row.names = FALSE)

# ---- compare ----
cols <- c("rwb_left", "rwb_right", "rwb_two", "ew_left", "ew_right", "ew_two", "b_ols", "b_ivx", "rho", "delta")
options(width = 200)
for (ret in ret_series) {
  r <- rep[rep$ret == ret, ]
  cmp <- merge(paper[, c("pred", cols)], r[, c("pred", cols)], by = "pred",
               suffixes = c("_paper", "_ivx"), sort = FALSE)
  cmp <- cmp[match(paper$pred, cmp$pred), ]
  for (cc in cols) cmp[[paste0(cc, "_diff")]] <- cmp[[paste0(cc, "_ivx")]] - cmp[[paste0(cc, "_paper")]]
  write.csv(cmp, file.path(here, paste0("table4-panelA-comparison-", ret, ".csv")), row.names = FALSE)

  cat("\n==================== return series:", ret, "====================\n")
  cat("\nAsymptotic Eicker-White p-values (paper | ivx)\n")
  print(cbind(cmp[, "pred", drop = FALSE], round(cmp[, grep("^ew_", names(cmp))], 3)), row.names = FALSE)
  cat("\nResidual wild bootstrap p-values (paper, B = 9999 | ivx, B =", B, ")\n")
  print(cbind(cmp[, "pred", drop = FALSE], round(cmp[, grep("^rwb_", names(cmp))], 3)), row.names = FALSE)
  cat("\nEstimates (paper | ivx; the paper scales tbl..infl by 100)\n")
  print(cbind(cmp[, "pred", drop = FALSE], round(cmp[, grep("^(b_|rho|delta)", names(cmp))], 3)), row.names = FALSE)
  cat("\nMean |diff|: EW p =", round(mean(abs(as.matrix(cmp[, grep("^ew_.*_diff$", names(cmp))]))), 3),
      "; RWB p =", round(mean(abs(as.matrix(cmp[, grep("^rwb_.*_diff$", names(cmp))]))), 3), "\n")
}

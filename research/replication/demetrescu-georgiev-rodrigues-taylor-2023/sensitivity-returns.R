# Sensitivity of the Table 4 (Panel A) asymptotic Eicker-White p-values to
# (i) the definition of the excess return and (ii) the Goyal-Welch data vintage.
# Fast (no bootstrap). Usage: Rscript sensitivity-returns.R

here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
pkg <- normalizePath(file.path(here, "..", "..", ".."))
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(pkg, quiet = TRUE) else library(ivx)

lag1 <- function(x) c(NA, head(x, -1))
build <- function(file) {
  gw <- read.csv(file.path(here, file), na.strings = c("NaN", ""), check.names = FALSE)
  gw$Index <- as.numeric(gsub(",", "", gw$Index))
  d <- data.frame(
    Date = as.Date(paste0(gw$yyyymm, "01"), "%Y%m%d"),
    vw_log    = log(1 + gw$CRSP_SPvw) - log(1 + gw$Rfree),    # paper's stated definition
    vw_simple = gw$CRSP_SPvw - gw$Rfree,
    vwx_log   = log(1 + gw$CRSP_SPvwx) - log(1 + gw$Rfree),   # ex-dividend
    idx_log   = log(gw$Index / lag1(gw$Index)) - log(1 + gw$Rfree),
    dp = log(gw$D12) - log(gw$Index), dy = log(gw$D12) - log(lag1(gw$Index)),
    ep = log(gw$E12) - log(gw$Index), de = log(gw$D12) - log(gw$E12),
    svar = gw$svar, bm = gw$`b/m`, ntis = gw$ntis, tbl = gw$tbl, lty = gw$lty, ltr = gw$ltr,
    tms = gw$lty - gw$tbl, dfy = gw$BAA - gw$AAA, dfr = gw$corpr - gw$ltr, infl = gw$infl
  )
  d[d$Date >= as.Date("1926-12-01") & d$Date <= as.Date("2020-12-01"), ]
}
paper <- c(dp = .5099, dy = .1492, ep = .1716, de = .6025, svar = .7712, bm = .5939, ntis = .2339,
           tbl = .0730, lty = .0543, ltr = .1639, tms = .5132, dfy = .9962, dfr = .3943, infl = .1113)
p_two <- function(d, ret, p) {
  unname(2 * (1 - pnorm(abs(ivx(as.formula(paste(ret, "~", p)), d, robust = TRUE)$tstat))))
}
rets <- c("vw_log", "vw_simple", "vwx_log", "idx_log")
out <- list()
for (v in c("goyal-welch-monthly.csv", "goyal-welch-monthly-2021.csv")) {
  d <- build(v)
  tab <- sapply(rets, function(r) sapply(names(paper), function(p) p_two(d, r, p)))
  out[[v]] <- cbind(paper = paper, round(tab, 3))
  cat("\n", v, " — two-sided Eicker-White p-values\n", sep = "")
  print(out[[v]])
  cat("mean |p - paper|:\n"); print(round(colMeans(abs(tab - paper)), 3))
}
write.csv(do.call(rbind, lapply(names(out), function(v) cbind(vintage = v, pred = rownames(out[[v]]), as.data.frame(out[[v]])))),
          file.path(here, "sensitivity-returns.csv"), row.names = FALSE)

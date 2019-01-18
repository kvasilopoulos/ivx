#' Unit root and stationarity test overview
#'
#'
#' @param y dataset
#' @param model include constant and or trend
#' @param lag.max maximum number of lags
#'
#' @importFrom urca ur.df ur.ers ur.pp ur.kpss
#' @importFrom dplyr select
#' @importFrom purrr reduce map
#' @importFrom stats lm.fit
#'
#' @export
#'
#' @examples
#' \donttest{
#' names <- datam %>%
#'   select(-Date, -Ret) %>%
#'   names()
#'
#' ur <- cbind(names, datam %>%
#'            select(-Date, -Ret) %>%
#'            map(uroot) %>%
#'            reduce(rbind))
#'}
uroot <- function(y, model = c("constant", "trend"), lag.max = 4) {

  nr <- length(y)
  yt <- as.matrix(y[2:nr])
  ylag <- cbind(1, y[1:(nr - 1)]) # 1 for constant
  Rn <- lm.fit(ylag, yt)$coefficients[[2]]

  # ADF Test
  adf <- ur.df(y, selectlags = c("BIC"))
  adf_tstat <- adf@teststat[1]
  adf_cv <- adf@cval

  # DF-GLS Test
  dfgls <- ur.ers(y)
  dfgls_tstat <- dfgls@teststat
  dfgls_cv <- dfgls@cval

  # PP-Test
  pp <- ur.pp(y, type = "Z-tau")
  pp_tstat <- pp@teststat
  pp_cv <- pp@cval

  # KPSS-Test
  kpss <- ur.kpss(y)
  kpss_tstat <- kpss@teststat
  kpss_cv <- kpss@cval

  output <- data.frame(Rn = Rn,
                       ADF = signify_star(adf_tstat, adf_cv),
                       DF_GLS = signify_star(dfgls_tstat, dfgls_cv),
                       PP = signify_star(pp_tstat, pp_cv),
                       KPSS = signify_star(kpss_tstat, kpss_cv))

  output
}

#' @importFrom dplyr case_when
signify_star <- function(tstat, cv, stat = TRUE) {
  tstat <- round(tstat, 3)
  if (stat) {
    nstars <- case_when(
      tstat < pct(cv, "10pct")  ~ 0,
      tstat >= pct(cv, "10pct") & tstat < pct(cv, "5pct") ~ 1,
      tstat >= pct(cv, "5pct") & tstat < pct(cv, "1pct") ~ 2,
      tstat >= pct(cv, "10pct") ~ 3
    )
  }else{
    # this is unit root significance starts
    nstars <- case_when(
      tstat < pct(cv, "10pct")  ~ 0,
      tstat >= pct(cv, "10pct") & tstat < pct(cv, "5pct") ~ 1,
      tstat >= pct(cv, "5pct") & tstat < pct(cv, "1pct") ~ 2,
      tstat >= pct(cv, "10pct") ~ 3
    )
  }

  if (nstars > 0) {
    output <-  paste0(tstat, paste(rep("*", nstars), collapse = ""))
  }else{
    output <- tstat
  }
 output
}

pct <- function(x, string) {
  x[which(colnames(x) == string)]
}

df_gls <- function(x) {
  nr <- NROW(x)
  rho <- 1 - 7.0/nr
  yt <- c(x[1], x[2:nr] - rho*x[1:(nr - 1)])
  xt <- c(1, rep(1 - rho, nr - 1))
  mu_reg <- summary(lm(yt ~ 0 + xt))
  yd <- yt - coef(mu_reg)[1]
  ydlag <-  yd[1:(nr - 1)]
  yddiff <- diff(yd)
  dfgls_reg <- summary(lm(yddiff ~ 0 + ydlag))
  tstat <- coef(dfgls_reg)[3]
  return(tstat)
}

#' Campel Yogo bonferroni confidence intervals Q-test
#'
#' @param y the regressand
#' @param x the regressor
#' @param method c("bonferroni", "infeasible", "sup")
#' @param cest for sup test the user have to provide an estimation of c
#' @param sig 10\%
#'
#' @importFrom stats rnorm
#' @return the upper and lower bound of the regression coefficient and the
#' autoregressive coefficient
#' @export
#'
#' @examples
#' x <- cumsum(rnorm(100))
#' y <- 2*x + rnorm(100)
#' qtest(y, x)
qtest <- function(y, x, method = c("bonferroni", "infeasible", "sup"),
                  cest = NULL, sig = 0.1) {

  pm <- pmatch(match.arg(method), c("bonferroni", "infeasible", "sup"))
  method <- c("bonferroni", "infeasible", "sup")[pm]

  # on.exit(return(method))

  x <- as.matrix(x)
  y <- as.matrix(y)
  nr <- NROW(x)

  xlag <- x[1:(nr - 1), , drop = FALSE]
  xt <- x[2:nr, , drop = FALSE]
  yt <- y[2:nr, , drop = F]

  nn <- NROW(xlag)
  l <- NCOL(xlag)

  main_reg <- summary(lm(yt ~ xlag))
  ar_reg <- summary(lm(xt ~ xlag))

  se_b <- coefficients(main_reg)[2,2]
  se_rho <- coefficients(ar_reg)[2,2]

  u <- residuals(main_reg)
  e <- residuals(ar_reg)
  nres <- NROW(u)
  sigma_u <- sqrt(crossprod(u)/(nres - 2))
  sigma_e <- sqrt(crossprod(e)/(nres - 2))
  sigma_ue <- t(u) %*% e/(nres - 2)
  # delta <- cor(u, e)
  # tstat and delta helps me to find the appropriate c_low and c_up
  delta <- sigma_ue/(sigma_u*sigma_e)
  tstat <- df_gls(x)


  # if (method == "ttest") {
  #   rho <- coefficients(ar_reg)[2,1]
  #   c_rho <- (rho - 1)*nn
  # }else if (method == "qtest") {
  #   c_low <- 1 # have to find a way to include them
  #   c_up <- 1
  # }else{
  #   rho <- 1
  # }

  if (method == "bonferroni") {

    delta_sim <- c(-0.999, seq(-0.975, -0.025, 0.025))
    dfgls_sim <- seq(1, -5, -0.1)
    # ts <- round(tstat, 1)

    s1 <- which.min(abs(dfgls_sim == tstat))
    s2 <- which.min(abs(delta_sim - delta[1]))

    c_low <- get("c_low_crit")[s1, s2]
    c_up <- get("c_up_crit")[s1, s2]

    rho_low <- 1 + c_low/nr
    rho_up <- 1 + c_up/nr
  } else if (method == "sup") {
    rho_low <- 1
    rho_up <- 1
  }else{
    if (is.null(cest)) {
      stop("Provide an estimation for c, in order to calculate the Q-test",
           call. = FALSE)
    }
    rho_low <- 1 + cest/nr
    rho_up <- 1 + cest/nr
  }

  zsig <- qnorm(1 - sig/2)

  # xlag_m <- xlag - mean(xlag)

  ylow <- yt - sigma_ue[1]*sigma_e[1]^(-2)*(xt - rho_low*xlag)
  main_low_reg <- coefficients(summary(lm(ylow ~ xlag))) # this should be xlag_m
  beta_low <- main_low_reg[2,1] - zsig*(1 - delta^2)^(1/2) * main_low_reg[2,2]

  yup <- yt - sigma_ue[1]*sigma_e[1]^(-2)*(xt - rho_up*xlag)
  main_up_reg <- coefficients(summary(lm(yup ~ xlag)))
  beta_up <- main_up_reg[2,1] + zsig*(1 - delta^2)^(1/2) * main_up_reg[2,2]

  return(
    structure(
      list(ci_b = c(beta_low, beta_up),
           ci_rho = c(rho_low, rho_up),
           ci_c = c(c_low, c_up),
           method = method,
           delta = delta[1],
           dfgls = tstat),
      class = "qtest")
    )
}

#' @export
print.qtest <- function(x, ...) {
  cat(" CI for beta: [", x$ci_b[1], ",", x$ci_b[2], "]")
  cat("\n CI for rho: [", x$ci_rho[1], ",", x$ci_rho[2], "]")
  cat("\n delta = ", x$delta)
}

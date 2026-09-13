#' Unified Empirical Likelihood Test for Predictability (Liu, Yang, Cai & Peng)
#'
#' `el_test` implements the unified empirical likelihood (EL) test of Liu,
#' Yang, Cai and Peng (2019, Section 2.2) for the predictive regression with
#' unknown intercept
#' \deqn{Y_t = \alpha + \beta_1 \Delta X_{t-1} + \beta_2 X_{t-2} + U_t,}
#' where the lagged difference of the predictor is included so that the
#' response can be stationary when the predictor is not. The intercept is
#' removed by differencing at lag \eqn{m = \lfloor n/2 \rfloor}
#' (\eqn{\tilde Y_t = Y_{t+m} - Y_t}, \eqn{\tilde X_t = X_{t+m} - X_t}), and the
#' EL function is built from the score equations
#' \eqn{\tilde Z_{t1} = \tilde e_t \Delta\tilde X_{t-1}} and
#' \eqn{\tilde Z_{t2} = \tilde e_t \tilde X_{t-2}/\sqrt{1 + \tilde X_{t-2}^2}},
#' \eqn{t = 3, \dots, m}, with \eqn{\tilde e_t} the model error. The weight on
#' the second equation makes its sample variance converge whether the predictor
#' is stationary, nearly integrated or a unit root, so that the profile EL
#' ratios for \eqn{H_0: \beta_2 = 0} (no predictability), \eqn{H_0: \beta_1 = 0}
#' and the joint null are \eqn{\chi^2(1)}, \eqn{\chi^2(1)} and \eqn{\chi^2(2)}
#' without knowing the persistence (Theorem 2). Single predictor only.
#'
#' @inheritParams ivx
#'
#' @return an object of class "el_test": a list with the EL ratio statistics
#' `stat` (named `beta2`, `beta1`, `joint`), their `p.value`, the
#' unconstrained EL estimates `estimate`, the OLS estimates `ols` of
#' \eqn{(\beta_1, \beta_2)} on the differenced data, and `m`.
#'
#' @references Liu, X., Yang, B., Cai, Z., & Peng, L. (2019). A unified test
#' for predictability of asset returns regardless of properties of predicting
#' variables. Journal of Econometrics, 208(1), 141-159.
#' @references Owen, A. B. (2001). Empirical Likelihood. Chapman & Hall.
#'
#' @export
#' @examples
#' el_test(Ret ~ DP, data = kms)
el_test <- function(formula, data, na.action) {
  cl <- match.call()
  fr <- ivx_frame(match.call(expand.dots = FALSE), parent.frame())
  if (NCOL(fr$x) != 1) stop("el_test() is defined for a single predictor", call. = FALSE)
  out <- el_test_fit(fr$y, drop(fr$x))
  out$call <- cl
  out$predictor <- colnames(fr$x)
  out
}

#' @rdname el_test
#' @param y response vector.
#' @param x predictor vector.
#' @export
el_test_fit <- function(y, x) {
  n <- length(x)
  if (length(y) != n) stop("incompatible dimensions")
  m <- floor(n / 2)
  if (m < 10) stop("too few observations", call. = FALSE)
  yt <- y[(1 + m):(2 * m)] - y[1:m]                 # tilde Y_t, t = 1..m
  xt <- x[(1 + m):(2 * m)] - x[1:m]
  t <- 3:m
  Y <- yt[t]
  DX <- xt[t - 1] - xt[t - 2]
  X2 <- xt[t - 2]
  W <- X2 / sqrt(1 + X2^2)
  scores <- function(b) {
    e <- Y - b[1] * DX - b[2] * X2
    cbind(e * DX, e * W)
  }
  ols <- lm.fit(cbind(DX, X2), Y)$coefficients
  names(ols) <- c("beta1", "beta2")

  el <- function(b) el_ratio(scores(b))
  # profiles: minimise over the nuisance coefficient from the OLS start
  p2 <- stats::optim(ols[1], function(b1) el(c(b1, 0)), method = "BFGS")
  p1 <- stats::optim(ols[2], function(b2) el(c(0, b2)), method = "BFGS")
  est <- stats::optim(ols, el, method = "BFGS")
  stat <- c(beta2 = p2$value, beta1 = p1$value, joint = el(c(0, 0)))
  pv <- 1 - stats::pchisq(stat, c(1, 1, 2))
  structure(
    list(stat = stat, p.value = pv, estimate = stats::setNames(est$par, c("beta1", "beta2")),
         ols = ols, m = m, n = n),
    class = "el_test"
  )
}

# -2 log empirical likelihood ratio that E[Z] = 0 (Owen, 2001, Ch. 3): Newton
# steps on the dual, with the quadratic-below-1/n pseudo-log so the objective
# is defined everywhere and large when 0 is outside the convex hull of the rows
el_ratio <- function(Z, maxit = 50, tol = 1e-9) {
  n <- nrow(Z)
  k <- ncol(Z)
  eps <- 1 / n
  logstar <- function(z) ifelse(z < eps, log(eps) - 1.5 + 2 * z / eps - z^2 / (2 * eps^2), log(pmax(z, eps)))
  d1 <- function(z) ifelse(z < eps, 2 / eps - z / eps^2, 1 / pmax(z, eps))
  d2 <- function(z) ifelse(z < eps, -1 / eps^2, -1 / pmax(z, eps)^2)
  lam <- rep(0, k)
  f <- function(l) -sum(logstar(1 + drop(Z %*% l)))
  val <- f(lam)
  for (i in seq_len(maxit)) {
    a <- 1 + drop(Z %*% lam)
    g <- -colSums(Z * d1(a))
    H <- crossprod(Z * sqrt(-d2(a)))
    step <- tryCatch(solve(H, g), error = function(e) g)
    s <- 1
    repeat {
      new <- lam - s * step
      nv <- f(new)
      if (nv <= val || s < 1e-8) break
      s <- s / 2
    }
    conv <- abs(val - nv) < tol
    lam <- new
    val <- nv
    if (conv) break
  }
  2 * sum(logstar(1 + drop(Z %*% lam)))
}

#' @rdname el_test
#' @param x an object of class "el_test".
#' @param digits minimal number of significant digits.
#' @param ... unused.
#' @export
print.el_test <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat("\nUnified empirical likelihood test (Liu, Yang, Cai & Peng, 2019), m = ", x$m, "\n\n", sep = "")
  tab <- cbind(Estimate = c(x$estimate, NA), "EL ratio" = x$stat[c("beta1", "beta2", "joint")],
               df = c(1, 1, 2), "Pr(> chi)" = x$p.value[c("beta1", "beta2", "joint")])
  rownames(tab) <- c("beta1 (dX[t-1])", "beta2 (X[t-2])", "joint")
  printCoefmat(tab, digits = digits, cs.ind = 1, tst.ind = 2:3, P.values = TRUE,
               has.Pvalue = TRUE, na.print = "", signif.stars = FALSE)
  cat("\n")
  invisible(x)
}

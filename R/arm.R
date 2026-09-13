#' Augmented Regression Method (Amihud, Hurvich & Wang)
#'
#' `arm` implements the multipredictor augmented regression method (mARM) of
#' Amihud, Hurvich and Wang (2009), a reduced-bias OLS alternative to IVX for
#' stationary but persistent predictors. A VAR(1) is fitted to the predictors,
#' its coefficient matrix is bias-corrected with the Nicholls and Pope (1988)
#' expansion (iterated), and the predictive regression is augmented with the
#' corrected VAR residuals, which removes the Stambaugh (1999) bias from the
#' slopes. Standard errors and the joint Wald statistic use the covariance
#' estimator of the paper (eqs 7-8), which adds the estimation uncertainty of
#' the VAR coefficients to the augmented-regression OLS variance.
#'
#' Unlike IVX the method assumes stationary predictors (all eigenvalues of the
#' VAR coefficient matrix inside the unit circle) and Gaussian innovations; it
#' is the natural benchmark for the "control function" approach of Elliott
#' (2011). Short horizon only.
#'
#' @inheritParams ivx
#' @param iter maximum number of bias-correction iterations (`K = 10` in the
#' paper); iteration stops earlier if the corrected VAR becomes non-stationary.
#'
#' @return an object of class `c("arm", "ivx")`, so the `ivx` methods apply.
#' Additional components: `phi` (coefficients on the augmentation residuals),
#' `Phi` (bias-corrected VAR(1) coefficient matrix, equations by row) and
#' `Phi_ols`.
#'
#' @references Amihud, Y., Hurvich, C. M., & Wang, Y. (2009).
#' Multiple-predictor regressions: Hypothesis testing. The Review of Financial
#' Studies, 22(1), 413-434.
#' @references Amihud, Y., & Hurvich, C. M. (2004). Predictive regressions: A
#' reduced-bias estimation method. Journal of Financial and Quantitative
#' Analysis, 39(4), 813-841.
#' @references Nicholls, D. F., & Pope, A. L. (1988). Bias in the estimation of
#' multivariate autoregressions. Australian Journal of Statistics, 30A, 296-309.
#'
#' @export
#' @examples
#' arm(Ret ~ DP, data = kms)
#'
#' summary(arm(Ret ~ DP + TBL, data = kms))
arm <- function(formula, data, iter = 10, na.action, contrasts = NULL,
                model = TRUE, x = FALSE, y = FALSE, ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  fr <- ivx_frame(match.call(expand.dots = FALSE), parent.frame(), contrasts)
  z <- arm_fit(fr$y, fr$x, iter = iter, ...)
  class(z) <- c("arm", "ivx")
  ivx_finish(z, fr, cl, model, ret.x, ret.y)
}

#' Fitter Function for the Augmented Regression Method
#'
#' Basic function called by `arm`. Should only be used directly by experienced
#' users.
#'
#' @inheritParams stats::lm.fit
#' @inheritParams arm
#' @param ... currently unused.
#' @export
#' @examples
#' arm_fit(kms$Ret, as.matrix(kms$DP))$tstat
arm_fit <- function(y, x, iter = 10, ...) {
  chkDots(...)
  n <- NROW(x)
  p <- NCOL(x)
  if (is.null(n)) stop("'x' must be a matrix")
  if (NROW(y) != n) stop("incompatible dimensions")
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("x", 1L:p)

  xt <- x[-1, , drop = FALSE]
  xl <- x[-n, , drop = FALSE]
  yt <- y[-1]
  nn <- n - 1

  # Step 1: OLS VAR(1) with intercept, Phi equations by row; Yule-Walker fallback
  X1 <- cbind(1, xl)
  XtXi <- solve(crossprod(X1))
  A <- XtXi %*% crossprod(X1, xt)
  Phi_ols <- t(A[-1, , drop = FALSE])
  if (!stationary(Phi_ols)) {
    xc <- scale(x, scale = FALSE)
    Phi_ols <- t(solve(crossprod(xc), crossprod(xc[-n, , drop = FALSE], xc[-1, , drop = FALSE])))
  }
  Sv <- function(Phi) {
    v <- xt - (rep(1, nn) %o% (colMeans(xt) - drop(Phi %*% colMeans(xl)))) - xl %*% t(Phi)
    crossprod(v) / nn
  }
  # Steps 2-3: iterate the Nicholls-Pope correction Phi_c = Phi_ols + b/n
  Phi <- Phi_ols
  S <- Sv(Phi)
  for (k in seq_len(iter)) {
    b <- np_bias(Phi, S)
    cand <- Phi_ols + b / nn
    if (!stationary(cand)) break
    Phi <- cand
    S <- Sv(Phi)
  }

  # Step 4: corrected residuals and the augmented regression
  mu_c <- colMeans(xt) - drop(Phi %*% colMeans(xl))
  vc <- xt - rep(1, nn) %o% mu_c - xl %*% t(Phi)
  Xa <- cbind(1, xl, vc)
  ols <- lm.fit(Xa, yt)
  cf <- ols$coefficients
  beta <- cf[2:(p + 1)]
  phi <- cf[(p + 2):(2 * p + 1)]
  s2e <- sum(ols$residuals^2) / (nn - (2 * p + 1))

  # eqs (7)-(8): cov[Phi_c' phi] from the OLS VAR covariance Sigma_v (x) (X'X)^-1,
  # plus E[B] from the partial residuals r_j of x_{j,t-1} on the other regressors
  V1 <- drop(crossprod(phi, S %*% phi)) * XtXi[-1, -1, drop = FALSE]
  R <- sapply(seq_len(p), function(j) lm.fit(Xa[, -(j + 1), drop = FALSE], xl[, j])$residuals)
  R <- matrix(R, ncol = p)
  d <- colSums(R^2)
  EB <- s2e * crossprod(R) / (d %o% d)
  V <- V1 + EB
  dimnames(V) <- list(cnames, cnames)

  se <- sqrt(diag(V))
  tstat <- beta / se
  names(beta) <- names(se) <- names(tstat) <- names(phi) <- cnames
  dimnames(Phi) <- dimnames(Phi_ols) <- list(cnames, cnames)
  fitted <- drop(xl %*% beta)
  u <- yt - cf[1] - fitted

  list(
    coefficients = beta,
    se = se,
    tstat = tstat,
    phi = phi,
    Phi = Phi,
    Phi_ols = Phi_ols,
    intercept = unname(cf[1]),
    fitted = fitted,
    residuals = u,
    Wald_Joint = drop(crossprod(beta, solve(V, beta))),
    Wald_Ind = tstat^2,
    rank = p,
    horizon = 1,
    df.residuals = nn - (2 * p + 1),
    df = p,
    assign = attr(x, "assign"),
    cnames = cnames,
    AR = data.frame(Rn = diag(Phi), row.names = cnames),
    delta = matrix(cor(u, vc), 1, p, dimnames = list(NULL, cnames)),
    vcov = V,
    robust = FALSE,
    tuning = list(iter = iter)
  )
}

stationary <- function(Phi) all(Mod(eigen(Phi, only.values = TRUE)$values) < 1)

# Nicholls & Pope (1988): E[Phi_hat - Phi] = -b/n with
# b = Sigma_v [ (I - Phi')^-1 + Phi'(I - Phi'^2)^-1 + sum_lambda lambda (I - lambda Phi')^-1 ] Gamma_0^-1,
# Gamma_0 = cov(x_t) implied by (Phi, Sigma_v) through the Lyapunov equation so that
# b is a function of the VAR parameters only (Section 2, step 3); reduces to
# 1 + 3 phi in the scalar case (Kendall, 1954)
np_bias <- function(Phi, S) {
  p <- nrow(Phi)
  I <- diag(p)
  Pt <- t(Phi)
  G0 <- matrix(solve(diag(p^2) - kronecker(Phi, Phi), as.vector(S)), p, p)
  G0i <- solve(G0)
  lam <- eigen(Phi, only.values = TRUE)$values
  spec <- Reduce(`+`, lapply(lam, function(l) l * solve(I - l * Pt)))
  b <- S %*% (solve(I - Pt) + Pt %*% solve(I - Pt %*% Pt) + spec) %*% G0i
  Re(b)
}

#' @rdname arm
#' @param x an object of class "arm".
#' @inheritParams stats::summary.lm
#' @export
print.arm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Augmented regression method (reduced-bias OLS)\n\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

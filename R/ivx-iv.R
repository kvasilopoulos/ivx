#' IV Predictability Tests with Combined Instruments
#'
#' `ivx_iv` implements the instrumental-variable tests of Breitung and
#' Demetrescu (2015): the predictive regression is estimated by 2SLS with
#' Eicker-White standard errors (their eq. 12) using instruments that are less
#' persistent than the predictor. Two families are available. Type-I
#' instruments are transformations of the predictor itself: the fractional
#' difference \eqn{\Delta_+^{d} x_{t-1}} (`"frac"`) and the long difference
#' \eqn{x_{t-1} - x_{t-1-k_T}} (`"diff"`), which keep power when the predictor
#' is stationary. Type-II instruments are deterministic and correlate with a
#' near-integrated predictor only: the sine function \eqn{\sin(\pi t/T)}
#' (`"sin"`). Their 2SLS combination (`"comb"`, the paper's `IVcomb` and the
#' authors' recommendation) is asymptotically dominated by whichever instrument
#' is informative, so the squared t-ratio (and the Wald statistic with several
#' predictors) is chi-square whatever the persistence of the predictor
#' (Theorems 3 and 7).
#'
#' With \eqn{K} predictors each type-I instrument is built per predictor and the
#' sine instruments use frequencies \eqn{\sin(k\pi t/T)}, \eqn{k = 1, \dots, K},
#' so that the instrument vector stays linearly independent (Assumption 5).
#' The IVX instrument of Kostakis et al. (2015) is the paper's "mild
#' integration" type-I case and is available through [ivx()].
#'
#' @inheritParams ivx
#' @param instruments instrument set; see Details.
#' @param d order of the fractional difference for `"frac"`, in (0, 1/2];
#' the paper uses 1/2.
#' @param kappa,eta the long-difference lag is \eqn{k_T = \lfloor \kappa T^\eta \rfloor}
#' (paper: 0.2 and 0.85), truncated to \eqn{t - 1}.
#'
#' @return an object of class `c("ivx_iv", "ivx")` so the `ivx` methods apply;
#' `instruments` holds the instrument matrix aligned with the regressors.
#'
#' @references Breitung, J., & Demetrescu, M. (2015). Instrumental variable and
#' variable addition based inference in predictive regressions. Journal of
#' Econometrics, 187(1), 358-375.
#'
#' @export
#' @examples
#' ivx_iv(Ret ~ DP, data = kms)
#' summary(ivx_iv(Ret ~ DP + TBL, data = kms, instruments = "frac"))
ivx_iv <- function(formula, data, instruments = c("comb", "sin", "frac", "diff"),
                   d = 0.5, kappa = 0.2, eta = 0.85, na.action, contrasts = NULL,
                   model = TRUE, x = FALSE, y = FALSE, ...) {
  instruments <- match.arg(instruments)
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  fr <- ivx_frame(match.call(expand.dots = FALSE), parent.frame(), contrasts)
  z <- ivx_iv_fit(fr$y, fr$x, instruments = instruments, d = d, kappa = kappa, eta = eta, ...)
  class(z) <- c("ivx_iv", "ivx")
  ivx_finish(z, fr, cl, model, ret.x, ret.y)
}

#' Fitter Function for IV Predictability Tests
#'
#' Basic function called by `ivx_iv`. Should only be used directly by
#' experienced users.
#'
#' @inheritParams stats::lm.fit
#' @inheritParams ivx_iv
#' @param ... currently unused.
#' @export
#' @examples
#' ivx_iv_fit(kms$Ret, as.matrix(kms$DP))$tstat
ivx_iv_fit <- function(y, x, instruments = "comb", d = 0.5, kappa = 0.2, eta = 0.85, ...) {
  chkDots(...)
  n <- NROW(x)
  l <- NCOL(x)
  if (is.null(n)) stop("'x' must be a matrix")
  if (NROW(y) != n) stop("incompatible dimensions")
  if (d <= 0 || d > 0.5) stop("`d` must be in (0, 1/2]", call. = FALSE)
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("x", 1L:l)

  xl <- x[-n, , drop = FALSE]                       # x_{t-1}, t = 2..n
  yt <- y[-1]
  nn <- n - 1

  Z <- switch(instruments,
    frac = frac_diff(xl, d),
    diff = long_diff(xl, kappa, eta),
    sin = sin_inst(nn, l),
    comb = cbind(sin_inst(nn, l), frac_diff(xl, d))
  )
  colnames(Z) <- switch(instruments,
    frac = paste0("frac_", cnames), diff = paste0("diff_", cnames),
    sin = paste0("sin", seq_len(l)), comb = c(paste0("sin", seq_len(l)), paste0("frac_", cnames))
  )

  # 2SLS with intercept in both stages (Section 3.2, f_t = 1)
  X1 <- cbind(1, xl)
  Z1 <- cbind(1, Z)
  ZZi <- solve(crossprod(Z1))
  P <- Z1 %*% ZZi %*% crossprod(Z1, X1)               # projected regressors
  A <- solve(crossprod(P, X1), crossprod(P, yt))
  coef <- drop(A)[-1]

  # Eicker-White sandwich of eq. (12) with OLS residuals of the predictive regression
  u <- lm.fit(X1, yt)$residuals
  Bi <- solve(crossprod(P, X1))
  V <- Bi %*% crossprod(P * u) %*% t(Bi)
  V <- V[-1, -1, drop = FALSE]
  dimnames(V) <- list(cnames, cnames)

  se <- sqrt(diag(V))
  tstat <- coef / se
  names(coef) <- names(se) <- names(tstat) <- cnames
  fitted <- drop(xl %*% coef)
  res <- yt - drop(A)[1] - fitted

  list(
    coefficients = coef,
    se = se,
    tstat = tstat,
    intercept = drop(A)[1],
    fitted = fitted,
    residuals = res,
    Wald_Joint = drop(crossprod(coef, solve(V, coef))),
    Wald_Ind = tstat^2,
    rank = l,
    horizon = 1,
    df.residuals = nn - l - 1,
    df = l,
    assign = attr(x, "assign"),
    cnames = cnames,
    instruments = Z,
    type = instruments,
    delta = matrix(cor(u, xl - rbind(0, xl[-nn, , drop = FALSE])), 1, l, dimnames = list(NULL, cnames)),
    vcov = V,
    robust = TRUE,
    tuning = list(d = d, kappa = kappa, eta = eta)
  )
}

# truncated fractional difference (1 - L)^d x_t = sum_{j=0}^{t-1} pi_j x_{t-j},
# pi_0 = 1, pi_j = pi_{j-1} (j - 1 - d) / j
frac_diff <- function(x, d) {
  n <- NROW(x)
  pi <- cumprod(c(1, (seq_len(n - 1) - 1 - d) / seq_len(n - 1)))
  z <- apply(x, 2, function(v) stats::convolve(v, rev(pi), type = "open")[seq_len(n)])
  matrix(z, ncol = NCOL(x))
}

# x_{t} - x_{t-k}, k = min(floor(kappa n^eta), t - 1); x_0 taken as 0 (so t = 1 gives x_1)
long_diff <- function(x, kappa, eta) {
  n <- NROW(x)
  k <- max(1L, floor(kappa * n^eta))
  lagk <- pmax(seq_len(n) - k, 0L)
  x0 <- rbind(0, x)
  x - x0[lagk + 1, , drop = FALSE]
}

sin_inst <- function(n, l) {
  t <- seq_len(n)
  sapply(seq_len(l), function(k) sin(k * pi * t / n))
}

#' @rdname ivx_iv
#' @param x an object of class "ivx_iv".
#' @inheritParams stats::summary.lm
#' @export
print.ivx_iv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("2SLS with instruments: ", paste(colnames(x$instruments), collapse = ", "), "\n\nCoefficients:\n", sep = "")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

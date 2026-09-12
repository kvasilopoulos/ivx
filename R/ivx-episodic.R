#' Subsample IVX Tests for Episodic Predictability
#'
#' Tests for "pockets" of predictability using the suprema of sequences of
#' subsample IVX statistics (Demetrescu et al., 2023, Section 3.2; Demetrescu
#' et al., 2022). For each window the IVX statistic is computed from the
#' window's observations with the full-sample instrument (eqs 15-17); the test
#' statistics are the maximum (right-tailed), minimum (left-tailed) and maximum
#' squared (two-sided) t-ratio over the sequence for a single predictor, and
#' the maximum Wald statistic for several predictors (Remark 11). P-values are
#' obtained by wild bootstrap (Algorithms 1-2), which the paper shows to be
#' asymptotically valid for these sup-functionals.
#'
#' @inheritParams ivx_boot
#' @param object an object of class "ivx" fitted with `horizon = 1`.
#' @param scheme `"rolling"` windows of fixed width, `"forward"` recursive
#' windows starting at the first observation, or `"backward"` recursive windows
#' ending at the last observation.
#' @param window fraction of the sample: the window width for `"rolling"`, the
#' warm-in fraction \eqn{\tau_L} for `"forward"`, and the latest start
#' \eqn{\tau_U} for `"backward"`.
#' @param robust logical; use Eicker-White standard errors in the subsample
#' statistics.
#' @param type bootstrap scheme, see [ivx_boot()]. The default fixed regressor
#' wild bootstrap is the scheme used by Demetrescu et al. (2022).
#'
#' @return an object of class "ivx_episodic": the observed sequence
#' (`sequence`, one row per window with its start/end index and statistics),
#' the sup statistics (`statistic`) and their bootstrap p-values (`p.value`),
#' plus the bootstrap draws (`boot`).
#'
#' @references Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
#' (2022). Testing for episodic predictability in stock returns. Journal of
#' Econometrics, 227(1), 85-113.
#' @references Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
#' (2023). Extensions to IVX methods of inference for return predictability.
#' Journal of Econometrics, 237(2), 105271.
#'
#' @export
#' @examples
#' mod <- ivx(Ret ~ DP, data = kms)
#' ivx_episodic(mod, scheme = "rolling", window = 0.2, B = 99, seed = 1)
ivx_episodic <- function(object, scheme = c("rolling", "forward", "backward"), window = 0.2,
                         robust = FALSE, B = 999, type = c("frwb", "rwb"), ar_max = 5,
                         dist = c("rademacher", "normal"), seed = NULL, cores = 1L) {
  scheme <- match.arg(scheme)
  type <- match.arg(type)
  dist <- match.arg(dist)
  if (object$horizon != 1) stop("subsample tests are available for horizon = 1 only", call. = FALSE)
  if (window <= 0 || window >= 1) stop("`window` must be in (0, 1)", call. = FALSE)
  if (!is.null(seed) && cores == 1L) set.seed(seed)

  bs <- boot_setup(object, type, ar_max, dist)
  n <- bs$n
  l <- bs$l
  nn <- n - 1
  rho <- 1 - object$tuning$cz / nn^object$tuning$beta

  # windows on the regression sample t = 2..n (0-based rows of y_t, x_{t-1}, z_{t-1})
  w <- floor(window * nn)
  if (w < l + 2) stop("`window` too small for the sample size", call. = FALSE)
  win <- switch(scheme,
    rolling  = cbind(start = 0:(nn - w), end = (w - 1):(nn - 1)),
    forward  = cbind(start = 0L, end = (w - 1):(nn - 1)),
    backward = cbind(start = 0:w, end = nn - 1)
  )

  stats <- function(y, x) {
    z <- ivx_instrument(x, rho)
    sub_ivx_cpp(y[-1], x[-n, , drop = FALSE], z, win[, "start"], win[, "end"], robust)
  }
  sup <- function(s) {
    if (l == 1) c(sup = max(s[, 2]), inf = min(s[, 2]), sup_sq = max(s[, 2]^2))
    else c(sup_wald = max(s[, 1]))
  }

  obs <- stats(bs$y, bs$x)
  stat <- sup(obs)

  run_chunk <- function(nb) {
    out <- matrix(NA_real_, nb, length(stat))
    for (b in seq_len(nb)) {
      d <- boot_sample(bs)
      out[b, ] <- sup(stats(d$y, d$x))
    }
    out
  }
  boot <- if (cores > 1L) boot_parallel(run_chunk, B, cores, seed) else run_chunk(B)
  colnames(boot) <- names(stat)
  p <- if (l == 1) {
    c(sup = mean(boot[, "sup"] >= stat["sup"]), inf = mean(boot[, "inf"] <= stat["inf"]),
      sup_sq = mean(boot[, "sup_sq"] >= stat["sup_sq"]))
  } else {
    c(sup_wald = mean(boot[, "sup_wald"] >= stat["sup_wald"]))
  }

  sequence <- data.frame(start = win[, "start"] + 2L, end = win[, "end"] + 2L, obs)
  names(sequence)[-(1:2)] <- c("Wald", paste0("t_", object$cnames))

  structure(
    list(
      call = object$call, scheme = scheme, window = window, robust = robust,
      type = type, B = B, cnames = object$cnames,
      sequence = sequence, statistic = stat, p.value = p, boot = boot
    ),
    class = "ivx_episodic"
  )
}

#' @rdname ivx_episodic
#' @param x an object of class "ivx_episodic".
#' @param digits minimal number of significant digits.
#' @param ... unused.
#' @export
print.ivx_episodic <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Subsample IVX tests, ", x$scheme, " scheme (window = ", x$window, "), ",
      nrow(x$sequence), " windows", if (x$robust) ", Eicker-White s.e." else "", "\n",
      switch(x$type, rwb = "Residual", frwb = "Fixed regressor"),
      " wild bootstrap, B = ", x$B, "\n\n", sep = "")
  lab <- if (length(x$cnames) == 1) {
    c(sup = "sup t   (H1: beta > 0)", inf = "inf t   (H1: beta < 0)", sup_sq = "sup t^2 (H1: beta != 0)")
  } else {
    c(sup_wald = "sup Wald")
  }
  tab <- data.frame(statistic = formatC(x$statistic, digits = digits),
                    "bootstrap p" = format.pval(x$p.value, digits = digits),
                    row.names = lab[names(x$statistic)], check.names = FALSE)
  print(tab)
  cat("\n")
  invisible(x)
}

# IVX instrument z_{t-1} for t = 2..n: z_1 = 0, z_t = rho z_{t-1} + dx_t
ivx_instrument <- function(x, rho) {
  n <- NROW(x)
  z <- apply(diff(x), 2, function(v) stats::filter(v, rho, "recursive"))
  rbind(0, matrix(z, ncol = NCOL(x))[-(n - 1), , drop = FALSE])
}

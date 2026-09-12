#' Wild Bootstrap Inference for IVX Models
#'
#' Computes bootstrap p-values for the IVX Wald and t statistics using the
#' residual wild bootstrap (RWB) or the fixed regressor wild bootstrap (FRWB) of
#' Demetrescu et al. (2023). The null hypothesis of no predictability is imposed
#' on the bootstrap samples. RWB rebuilds the regressor from an AR fit and its
#' residuals (multiplied by the same wild multiplier as the predictive-regression
#' residuals, so the innovation correlation is preserved); FRWB keeps the
#' regressors and instruments fixed and only resamples the response.
#'
#' @param object an object of class "ivx" (not "ivx_ar"), fitted without weights.
#' @param B number of bootstrap replications.
#' @param type bootstrap scheme: `"rwb"` (residual wild bootstrap, recommended
#' for strongly persistent regressors) or `"frwb"` (fixed regressor wild bootstrap).
#' @param ar_max maximum lag order of the (vector) autoregression fitted to the
#' regressors by the RWB scheme; the order is selected by BIC (Remark 24).
#' @param dist distribution of the wild multipliers.
#' @param seed optional integer seed. With `cores > 1` the L'Ecuyer-CMRG
#' streams of the \pkg{parallel} package are used, so results are reproducible
#' for a given `seed` and `cores` but differ from the serial run.
#' @param cores number of CPU cores. Uses forking on Unix and a PSOCK cluster on
#' Windows (the package must be installed for the workers to load it).
#'
#' @return an object of class "ivx_boot": a list with the observed statistics,
#' the bootstrap distributions (`boot`), and bootstrap p-values (`p.value`).
#' `p.value$tstat` has one column per alternative: two-sided, `beta < 0` and
#' `beta > 0`.
#'
#' @references Demetrescu, M., Georgiev, I., Rodrigues, P. M. M., & Taylor, A. M. R.
#' (2023). Extensions to IVX methods of inference for return predictability.
#' Journal of Econometrics, 237(2), 105271.
#'
#' @export
#' @importFrom stats embed lm.fit model.frame model.matrix model.response rnorm
#' @examples
#' mod <- ivx(Ret ~ DP + TBL, data = kms)
#' ivx_boot(mod, B = 199, seed = 1)
ivx_boot <- function(
  object,
  B = 999,
  type = c("rwb", "frwb"),
  ar_max = 5,
  dist = c("rademacher", "normal"),
  seed = NULL,
  cores = 1L
) {
  type <- match.arg(type)
  dist <- match.arg(dist)
  if (!inherits(object, "ivx") || inherits(object, "ivx_ar")) {
    stop(
      "`object` must be of class 'ivx' (ivx_ar is not supported)",
      call. = FALSE
    )
  }
  if (!is.null(object$weights)) {
    stop("bootstrap is not available for weighted fits", call. = FALSE)
  }
  if (!is.null(seed) && cores == 1L) {
    set.seed(seed)
  }

  x <- model.matrix(object)
  y <- model.response(model.frame(object), "numeric")
  if (!is.null(object$offset)) {
    y <- y - object$offset
  }
  n <- NROW(x)
  l <- NCOL(x)
  tun <- object$tuning
  fit <- function(y, x) {
    ivx_fit_cpp(
      y,
      x,
      object$horizon,
      tun$beta,
      tun$cz,
      tun$bandwidth,
      object$robust
    )
  }

  # Step 1: predictive regression residuals, t = 2, ..., n
  u_hat <- object$ols$residuals

  # Step 2 (RWB): VAR(q) on the regressors (Remark 22), residuals aligned with
  # u_hat (t = 2, ..., n), zero for t <= q
  if (type == "rwb") {
    va <- var_ols(x, ar_max)
    v_hat <- matrix(0, n - 1, l)
    v_hat[va$q:(n - 1), ] <- va$resid
  }

  draw <- if (dist == "rademacher") {
    function(m) sample(c(-1, 1), m, replace = TRUE)
  } else {
    function(m) rnorm(m)
  }

  # one chunk of replications; returns B_k x (1 + 2 l) matrix of statistics
  run_chunk <- function(nb) {
    out <- matrix(NA_real_, nb, 1 + 2 * l)
    for (b in seq_len(nb)) {
      r <- draw(n - 1)
      # Steps 3-4: impose the null, y* = u*; rebuild x* (RWB) or keep it fixed (FRWB)
      y_b <- c(0, r * u_hat)
      x_b <- if (type == "rwb") var_sim_cpp(rbind(0, r * v_hat), va$A) else x
      z <- fit(y_b, x_b)
      out[b, ] <- c(z$wivx, z$wivxind, z$zinvxind)
    }
    out
  }
  stats <- if (cores > 1L) boot_parallel(run_chunk, B, cores, seed) else run_chunk(B)
  boot_joint <- stats[, 1]
  boot_ind <- stats[, 1 + seq_len(l), drop = FALSE]
  boot_t <- stats[, 1 + l + seq_len(l), drop = FALSE]
  colnames(boot_ind) <- colnames(boot_t) <- object$cnames

  # Step 6: bootstrap p-values
  p_t <- cbind(
    "two-sided" = colMeans(abs(boot_t) >= rep(abs(object$tstat), each = B)),
    "less" = colMeans(boot_t <= rep(object$tstat, each = B)),
    "greater" = colMeans(boot_t >= rep(object$tstat, each = B))
  )
  rownames(p_t) <- object$cnames

  structure(
    list(
      call = object$call,
      type = type,
      B = B,
      dist = dist,
      coefficients = object$coefficients,
      tstat = object$tstat,
      Wald_Ind = object$Wald_Ind,
      Wald_Joint = object$Wald_Joint,
      boot = list(Wald_Joint = boot_joint, Wald_Ind = boot_ind, tstat = boot_t),
      p.value = list(
        Wald_Joint = mean(boot_joint >= object$Wald_Joint),
        Wald_Ind = colMeans(boot_ind >= rep(object$Wald_Ind, each = B)),
        tstat = p_t
      )
    ),
    class = "ivx_boot"
  )
}

#' @rdname ivx_boot
#' @param x an object of class "ivx_boot".
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to [printCoefmat()].
#' @export
print.ivx_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(
    switch(x$type, rwb = "Residual", frwb = "Fixed regressor"),
    " wild bootstrap, B = ",
    x$B,
    "\n\n",
    sep = ""
  )
  tab <- cbind(
    Estimate = x$coefficients,
    "t value" = x$tstat,
    "Wald Ind" = x$Wald_Ind,
    "Pr(> chi)" = x$p.value$Wald_Ind,
    "Pr(t < 0)" = x$p.value$tstat[, "less"],
    "Pr(t > 0)" = x$p.value$tstat[, "greater"]
  )
  cat("Coefficients (bootstrap p-values):\n")
  printCoefmat(
    tab,
    digits = digits,
    cs.ind = 1,
    tst.ind = 2:3,
    P.values = TRUE,
    has.Pvalue = TRUE,
    signif.stars = FALSE,
    ...
  )
  cat(
    "\nJoint Wald statistic: ",
    formatC(x$Wald_Joint, digits = digits),
    ", bootstrap p-value ",
    format.pval(x$p.value$Wald_Joint, digits = digits),
    "\n\n",
    sep = ""
  )
  invisible(x)
}

# OLS VAR(q) in levels (with intercept unless intercept = FALSE). `lag` is either the maximum order, with q
# chosen by `ic` on a common sample, or the fixed order when `fixed = TRUE`.
# Returns lag matrices A[[j]] (equation by row), q, and residuals for t = q+1..n.
var_ols <- function(x, lag, ic = c("bic", "aic"), fixed = FALSE, intercept = TRUE) {
  ic <- match.arg(ic)
  l <- NCOL(x)
  X <- function(e, q) {
    lags <- e[, l + seq_len(q * l), drop = FALSE]
    if (intercept) cbind(1, lags) else lags
  }
  if (fixed) {
    q <- lag
  } else {
    e <- embed(x, lag + 1)
    m <- nrow(e)
    pen <- if (ic == "bic") log(m) else 2
    crit <- sapply(seq_len(lag), function(q) {
      res <- lm.fit(X(e, q), e[, seq_len(l)])$residuals
      log(det(crossprod(res) / m)) + (intercept + q * l) * l * pen / m
    })
    q <- which.min(crit)
  }
  e <- embed(x, q + 1)
  fit <- lm.fit(X(e, q), e[, seq_len(l)])
  coef <- as.matrix(fit$coefficients)
  A <- lapply(seq_len(q), function(j) {
    t(coef[intercept + (j - 1) * l + seq_len(l), , drop = FALSE])
  })
  list(q = q, A = A, resid = as.matrix(fit$residuals))
}

# run `run_chunk` over B replications split across `cores` workers
boot_parallel <- function(run_chunk, B, cores, seed) {
  chunks <- tabulate(cut(seq_len(B), cores, labels = FALSE), cores)
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, loadNamespace("ivx"))
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    res <- parallel::parLapply(cl, chunks, run_chunk)
  } else {
    if (!is.null(seed)) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
    }
    res <- parallel::mclapply(chunks, run_chunk, mc.cores = cores, mc.set.seed = TRUE)
  }
  do.call(rbind, res)
}

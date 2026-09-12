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
#' @param ar_max maximum AR order (selected by AIC) for the regressor
#' autoregression used by the RWB scheme.
#' @param dist distribution of the wild multipliers.
#' @param seed optional integer passed to [set.seed()].
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
#' @importFrom stats ar.ols filter model.frame model.matrix model.response rnorm
#' @examples
#' mod <- ivx(Ret ~ DP + TBL, data = kms)
#' ivx_boot(mod, B = 199, seed = 1)
ivx_boot <- function(object, B = 999, type = c("rwb", "frwb"), ar_max = 5,
                     dist = c("rademacher", "normal"), seed = NULL) {
  type <- match.arg(type)
  dist <- match.arg(dist)
  if (!inherits(object, "ivx") || inherits(object, "ivx_ar")) {
    stop("`object` must be of class 'ivx' (ivx_ar is not supported)", call. = FALSE)
  }
  if (!is.null(object$weights)) {
    stop("bootstrap is not available for weighted fits", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  x <- model.matrix(object)
  y <- model.response(model.frame(object), "numeric")
  if (!is.null(object$offset)) y <- y - object$offset
  n <- NROW(x)
  l <- NCOL(x)
  tun <- object$tuning
  fit <- function(y, x) {
    ivx_fit_cpp(y, x, object$horizon, tun$beta, tun$cz, tun$bandwidth, object$robust)
  }

  # Step 1: predictive regression residuals, t = 2, ..., n
  u_hat <- object$ols$residuals

  # Step 2 (RWB): AR(p) on each regressor, residuals aligned with u_hat, NA -> 0
  if (type == "rwb") {
    ar_fit <- lapply(seq_len(l), function(j) {
      f <- ar.ols(x[, j], aic = TRUE, order.max = ar_max)
      if (f$order == 0) f <- ar.ols(x[, j], aic = FALSE, order.max = 1)
      f
    })
    ar_coef <- lapply(ar_fit, function(f) as.numeric(f$ar))
    v_hat <- sapply(ar_fit, function(f) {
      r <- as.numeric(f$resid)[-1]
      r[is.na(r)] <- 0
      r
    })
    v_hat <- matrix(v_hat, ncol = l)
  }

  draw <- if (dist == "rademacher") {
    function(m) sample(c(-1, 1), m, replace = TRUE)
  } else {
    function(m) rnorm(m)
  }

  boot_joint <- numeric(B)
  boot_ind <- matrix(NA_real_, B, l)
  boot_t <- matrix(NA_real_, B, l)
  for (b in seq_len(B)) {
    r <- draw(n - 1)
    # Steps 3-4: impose the null, y* = u*; rebuild x* (RWB) or keep it fixed (FRWB)
    y_b <- c(0, r * u_hat)
    x_b <- if (type == "rwb") {
      sapply(seq_len(l), function(j) {
        as.numeric(filter(c(0, r * v_hat[, j]), ar_coef[[j]], method = "recursive"))
      })
    } else {
      x
    }
    z <- fit(y_b, matrix(x_b, ncol = l))
    boot_joint[b] <- drop(z$wivx)
    boot_ind[b, ] <- drop(z$wivxind)
    boot_t[b, ] <- drop(z$zinvxind)
  }
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
      type = type, B = B, dist = dist,
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
  cat(switch(x$type, rwb = "Residual", frwb = "Fixed regressor"),
      " wild bootstrap, B = ", x$B, "\n\n", sep = "")
  tab <- cbind(
    Estimate = x$coefficients,
    "t value" = x$tstat,
    "Wald Ind" = x$Wald_Ind,
    "Pr(> chi)" = x$p.value$Wald_Ind,
    "Pr(t < 0)" = x$p.value$tstat[, "less"],
    "Pr(t > 0)" = x$p.value$tstat[, "greater"]
  )
  cat("Coefficients (bootstrap p-values):\n")
  printCoefmat(tab, digits = digits, cs.ind = 1, tst.ind = 2:3, P.values = TRUE,
               has.Pvalue = TRUE, signif.stars = FALSE, ...)
  cat("\nJoint Wald statistic: ", formatC(x$Wald_Joint, digits = digits),
      ", bootstrap p-value ", format.pval(x$p.value$Wald_Joint, digits = digits),
      "\n\n", sep = "")
  invisible(x)
}

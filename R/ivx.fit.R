ivx.fit <- function(y, x, h = 1, offset = NULL, ...) {

  if (is.null(n <- NROW(x))) stop("'x' must be a matrix")
  if (n == 0L)  stop("0 (non-NA) cases")
  p <- NCOL(x)
  if (p == 0L) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, rank = 0, df.residual = length(y)))
  }
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1)
    y <- drop(y)
  if (!is.null(offset))
    y <- y - offset
  if (NROW(y) != n)  stop("incompatible dimensions")

  chkDots(...)

  z <- ivx_fit_cpp(y, x, h)
  coef <- drop(z$Aivx)
  coef_ols <- t(z$Aols)

  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("x", 1L:p)
  # nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
  z$coefficients <- coef
  # r1 <- y - z$residuals
  if (!is.null(offset)) r1 <- r1 + offset
  if (is.matrix(y)) {
    dimnames(coef) <- list(cnames, colnames(y))
    # dimnames(z$effects) <- list(nmeffects, colnames(y))
  }else{
    names(coef) <- cnames
  }

  output <- structure(list(coefficients =  coef,
                           Wald_Joint = z$wivx,
                           Wald_Ind = z$wivxind,
                           horizon = h,
                           df = z$df,
                           cnames = cnames,
                           AR = data.frame(Rn = z$Rn,
                                           Rz = z$Rz,
                                           row.names = cnames),
                           delta = z$delta,
                           varcov = z$varcov,
                           coefficients_ols = coef_ols
                           ))


  output
}






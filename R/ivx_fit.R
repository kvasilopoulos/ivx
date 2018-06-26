ivx.fit <- function(y, x, h = 1) {

  x <- as.matrix(x)
  y <- as.matrix(y)
  nr <- NROW(x)

  xlag <- x[1:(nr - 1), , drop = FALSE]
  xt <- x[2:nr, , drop = FALSE]
  yt <- y[2:nr, , drop = F]

  nn <- NROW(xlag)
  l <- NCOL(xlag)

  Wivx <- matrix(0, 2, 1)
  WivxInd <- matrix(0, 2, l)

  # predictive regression residual estimation
  lm1 <- lm(yt ~ xlag)
  Aols <- coefficients(lm1)
  epshat <- matrix(residuals(lm1))

  rn <- matrix(0, l, l)
  for (i in 1:l) {
    rn[i, i] <- lm(xt[, i] ~ 0 + xlag[, i])$coefficients
  }

  # autoregressive residual estimation
  u <- xt - xlag %*% rn

  # residuals' correlation matrix
  # corrmat <- cor(cbind(epshat, u))
  corrmat <- c(epshat, u)

  # covariance matrix estimation (predictive regression)
  covepshat <- crossprod(epshat) / nn

  # %covariance matrix estimation (autoregression)
  covu <- matrix(0, l, l)
  for (i in 1:nn) {
    covu <- covu + crossprod(u[i, , drop = FALSE])
  }
  covu <- covu / nn

  # %covariance matrix between 'epshat' and 'u'
  covuhat <- matrix(0, 1, l)
  for (i in 1:l) {
    covuhat[, i] <- sum(epshat * u[, i])
  }
  covuhat <- t(covuhat) / nn


  m <- floor(nn^(1 / 3))
  uu <- matrix(0, l, l)
  for (i in 1:m) {
    a <- matrix(0, l, l)
    for (j in (i + 1):nn) {
      a <- a + t(u[j, , drop = F]) %*% u[j - i, , drop = F]
    }
    uu <- uu + (1 - i / (m + 1)) * a
  }
  uu <- uu / nn
  Omegauu <- covu + uu + t(uu)

  q <- matrix(0, m, l)
  for (i in 1:m) {
    p <- matrix(0, nn - i, l)
    for (j in (i + 1):nn) {
      p[j - i, ] <- u[j, , drop = F] * epshat[j - i]
    }
    q[i, ] <- (1 - i / (1 + m)) * colSums(p)
  }
  residue <- apply(q, 2, sum) / nn
  Omegaeu <- covuhat + residue

  # instrument construction
  n <- nn - h + 1
  Rz <- (1 - 1 / (nn^0.95)) * diag(l)
  diffx <- xt - xlag
  z <- matrix(0, nn, l)
  z[1, ] <- diffx[1, ]

  for (i in 2:nn) {
    z[i, ] <- z[i - 1, ] %*% Rz + diffx[i, ]
  }
  Z <- rbind(
    matrix(0, 1, l),
    z[1:(n - 1), , drop = F]
  )

  zz <- rbind(
    matrix(0, 1, l),
    z[1:(nn - 1), , drop = F]
  )


  ZK <- matrix(0, n, l)
  for (i in 1:(n - h + 1)) {
    ZK[i, ] <- colSums(zz[i:(i + h - 1), , drop = F])
  }

  yy <- matrix(0, n, 1)
  for (i in 1:(n - h + 1)) {
    yy[i] <- sum(yt[i:(i + h - 1), drop = F])
  }

  xK <- matrix(0, n, l)
  for (i in 1:(n - h + 1)) {
    xK[i, ] <- colSums(xlag[i:(i + h - 1), , drop = F])
  }

  meanxK <- colMeans(xK)
  Yt <- yy - mean(yy)
  Xt <- matrix(0, n, l)
  for (i in 1:l) {
    Xt[, i] <- xK[, i, drop = F] - meanxK[i] * matrix(1, n, 1)
  }

  library(pracma)
  # library(MASS)
  # tol <- max(size(A)) * eps(norm(A))

  Aivx <- t(Yt) %*% Z %*% pinv(t(Xt) %*% Z)
  meanzK <- colMeans(ZK)

  FM <- covepshat - t(Omegaeu) %*% Omegauu^(-1) %*% Omegaeu
  # meanzK drop to column dimension so I have to suse tcross instead of cross
  M <- crossprod(ZK) * covepshat[1] - n * tcrossprod(meanzK) * FM[1]

  H <- matrix(1, l, l)
  Q <- H %*% pinv(t(Z) %*% Xt) %*% M %*% pinv(t(Xt) %*% Z) * t(H)

  Wivx[1, 1] <- t(H %*% t(Aivx)) %*% pinv(Q) %*% (H %*% t(Aivx))
  Wivx[2, 1] <- 1 - pchisq(Wivx[1, 1], l)
  WivxInd[1, ] <- Aivx / t(diag(Q)^(1 / 2))
  WivxInd[2, ] <- 1 - pchisq(WivxInd[1, ]^2, 1)

  return(list(coefficients =  Aivx[-1], Wivx = Wivx,
              WivxInd = WivxInd, q = Q, corrmat))
}

print.ivx <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2L, quote = FALSE)
  } else {
    cat ("No coefficients\n")
  }
  invisible(x)
}



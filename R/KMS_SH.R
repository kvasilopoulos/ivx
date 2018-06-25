

lcov_mat1 <- function(x, lam){
	x <- as.matrix(x)
	nr <- NROW(x)
	nc <- NCOL(x)
	#lam <- floor(4*(nr/100)^(2/9)) # truncation parameter, see Newey West 1987
	for (i in 0:lam) {
		gamma <- (1 - i/(lam + 1))/nr*t(x[1:(nr - i), ]) %*% x[(1 + i):nr, ] # Bartlett window
		if (i == 0) {
			long <- gamma # Omega
			onesided <- gamma # Delta in KMS notation
			lambda <- 0
			sigma <- gamma
		}else {
			long <- long + gamma + t(gamma)
			onesided <- onesided + gamma
			lambda <- lambda + gamma
		}
	}
	return(list(long = long, onesided = onesided, lambda = lambda, sigma = sigma))
}



kms <- function(yt, xt, constant){

  n1 <- NROW(xt)
  beta <- 0.95
  npow <- 1/3
  cz <- -1

  ydata <- yt[2:n1]
  xdata <- xt[1:(n1 - 1)]
  xdata2 <- xt[2:n1]
  dxdata <- xdata2 - xdata

  n <- length(xdata)
  p <- 2
  Rnhat <- 0
  RnhatCterm <- 0
  zinstr <- rep(0,n)

  ## OLS estimation of the 1st equation with constant: y{t} <- x{t-1} ##

  olsreg1 <- lm(ydata~xdata)

  u0 <- as.vector(olsreg1$resid)
  Alpha1C <- olsreg1$coef
  Alpha1 <- Alpha1C[2]

  ## OLS estimation of the 2nd equation: x{t} <- C*x{t-1}, C being diagonal ##
  ## Perhaps Change to include an intercept

  olsreg2 <- lm(xdata2~xdata-1)
  Rnhat <- olsreg2$coef[1]
  #print(Rnhat)
  # RnhatCterm <- olsreg2$coef[1]
  ux <- as.vector(olsreg2$resid)


  ## var-cov  estimation ##
  ut <- cbind(u0,ux)
  M1 <- n^npow
  M <- trunc(M1)
  vcovmat <- lcov_mat1(x = ut, lam = M) # returns long, onesided, sigma, and lambda

  lambda <- vcovmat$lambda
  ## ok
  lambdax0 <- lambda[1,2] #mine is the transpose of MKS [(ky+1):p,1:ky]

  sigmaut <- vcovmat$sigma
  omegamatrix <- vcovmat$long

  omegaxx <- omegamatrix[2,2]
  omega00b <- sigmaut[1,1]
  omega0xb <- sigmaut[1,2]
  omegafmd <- omega00b - (omega0xb + lambdax0) * (omega0xb + lambdax0)/omegaxx

  ## loop for the creation of the IVX instruments ##

  Rnz <- 1 + cz/(n^beta)
  zinstr1 <- stats::filter(dxdata,Rnz,"rec")
  zinstr4 <- rep(0,n)
  zinstr4[2:n] <- zinstr1[1:(n - 1)]

  zl1 <- zinstr4

  ## A~ matrix calculation ##

  ydatadm <- ydata - mean(ydata)
  xdatadm <- xdata - mean(xdata)

  Alphatilb <- sum(ydatadm*zl1)/sum(xdatadm*zl1)

  ##

  zl1bar <- mean(zl1)
  #Pz <- zl1(sum(zl1*zl1)*zl1'

  valphatilb <- Alphatilb#vec(Alphatilb)

  ## Wald Statistic

  MM <- sum(zl1*zl1)*omega00b - n*zl1bar*zl1bar*omegafmd
  QH <- MM/(sum(zl1*xdatadm)^2)
  Wn3 <-  valphatilb*valphatilb/QH
  pval <- 1 - pchisq(Wn3,1)
  return(list(Aivx = Alphatilb,WaldIVX = Wn3,pval = pval))
}


#,u0 <- u0,ux <- ux,omega0xb <- omega0xb,lambdax0 <- lambdax0,omegafmd <- omegafmd,zl1 <- zl1







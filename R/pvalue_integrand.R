integrand <- function(rtilde_b, ttilde_gg, r_g, r_bb, r_gg, z_bb, rhobarhat, k) {

  tau <- sqrt(-log(1 - max(ttilde_gg)))
  A <- B <- matrix(NA, 2, 2)

  denom <- cosh(2*tau) + cos(2*tau)
  num <- sinh(2*tau) + sin(2*tau)

  A[1, 1] <- (1 / tau) * num / denom
  A[1, 2] <- (1 / (tau^2)) * (2*sinh(tau)*sin(tau)) / denom
  A[2, 1] <- A[1, 2]
  A[2, 2] <- (1 / (2*(tau^3))) * num / denom

  B[1, 1] <- (1 / tau) * num / denom
  B[1, 2] <- (1 / (tau^2)) * (1 - (2*cosh(tau)*cos(tau)) / denom)
  B[2, 1] <- B[1,2]
  B[2, 2] <- (1 / (tau^2)) * (1 - (1 / (2*tau)) * num / denom)

  B_inv <- solve(B)

  remat <- -pracma::inv(A %*% B_inv %*% A + B) %*% A %*% B_inv
  imat <- B_inv + B_inv %*% A %*% remat

  remat11   <- remat[1, 1]
  remat1221 <- remat[1, 2] + remat[2, 1]
  remat22   <- remat[2, 2]

  imat11   <- imat[1, 1]
  imat1221 <- imat[1, 2] + imat[2, 1]
  imat22   <- imat[2, 2]

  # Multipliers
  multi <- 1 / (sign(cos(tau/2)) * sqrt(cosh(sqrt(-2*i*(tau^2)))) * sqrt(det(A + i*B)))
  multi_re <- Re(multi)
  multi_im <- Im(multi)

  r_b <- ((-log(1 - rtilde_b))^k - 2*r_g - 1)/(2*rhobarhat)
  q_g <- (-log(1 - rtilde_b))^k


# minus -------------------------------------------------------------------

  # Exponent minus real
  exp_m_real <- -(r_b^2)/(2*r_bb) + remat_11*q_g - remat_1221*sqrt(q_g)*z_bb + remat_22*(z_bb^2)
  # Exponent minus imaginary
  exp_m_imag <- log(1 - ttilde_gg)*r_gg + imat_11*q_g - imat_1221*sqrt(q_g)*z_bb + imat_22*(z_bb^2)

  f_minus <- exp(exp_m_real) * (multi_real * cos(exp_m_imag) - multi_imag * sin(exp_m_im)) * ((-log(1 - rtilde_b))^(k/2 - 1) / ((1 - rtilde_b)*(1 - ttilde_gg))) / (abs(rhobarhat)*sqrt(r_bb)*z_bb)


# plus --------------------------------------------------------------------

  exp_p_re <- -(r_b^2)/(2*r_bb) + remat_11*q_g + remat_1221*sqrt(q_g)*z_bb + remat_22*(z_bb^2)

  exp_p_im <-  log(1 - ttilde_gg)*r_gg + imat_11*q_g + imat_1221*sqrt(q_g)*z_bb + imat_22*(z_bb^2)

  f_plus <- exp(exp_p_re) * (multi_re * cos(exp_p_im) - multi_im * sin(exp_p_im)) * ((-log(1 - rtilde_b))^(k/2 - 1) / ((1 - rtilde_b)*(1 - ttilde_gg))) / (abs(rhobarhat)*sqrt(r_bb)*z_bb)

  return(list(minus = f_minus,
              plus = f_plus))

}


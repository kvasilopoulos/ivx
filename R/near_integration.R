
# % Parameterization of Monte Carlos
mc <- function(n = 500, nrep = 1000){
  eta = 0.05
  beta_0 = 0
  b_grid = c(0,5, 10, 15, 20)
  c_grid = c(0, -5, -10, -15, -20)
  rho_grid = c(-.5,.5)


  # % Configuration of procedure
  k = 4
  Itegral_upper_truncation_rtilde = exp(-(4.^(2./k))) #% Upper boundary in rtilde units is 1 - Integral_upper_truncation_rtilde
  Integral_lower_truncation_rtilde = 0.e-6 #% Lower boundary in rtilde units is Integral_lower_truncation_rtilde
  Integral_upper_truncation_ttilde = 1.e-5 #% Upper boundary in ttilde units is 1 - Integral_upper_truncation_ttilde
  Integral_lower_truncation_ttilde = 0.e-6 #% Lower boundary in ttilde units is Integral_lower_truncation_ttilde

  # Initialization

  pvalue <- Integral_full <- Integral_full_plus  <- Integral_full_minus <-
    Integral_lower <- Integral_lower_plus <- Integral_lower_minus <-
    Integral_upper <- Integral_upper_plus <- Integral_upper_minus <-
    Integral_reliability <- Integral_reliability_plus <-
    Integral_reliability_minus <- rtilde_b_thresholds <-
    matrix(0, nrow = n, ncol = length(b_grid))

  x = matrix(0, nrep, 1)
  constant = matrix(1, nrep ,1)
  for (rho_index in 1:length(rho_grid)) {
    rho  <-  rho_grid[rho_index]

    for (c_index in 1:length(c_grid)) {
      c = c_grid[c_index]
      gamma = 1 + c/nrep # True value of gamma

      for (b_index in 1:length(b_grid)) {
        beta = beta_0 + sqrt(1 - rho^2)*b_grid[b_index]/T # True value of beta

        for (i in 1:n) { # Generate data set
          eps_x = rnorm(nrep, 1)
          eps_yy_x = sqrt(1 - rho^2)*rnorm(nrep, 1)
          eps_y = rho*eps_x + eps_yy_x

          x[1] = eps_x[1]
          for (t in 2:nrep) x[t] = gamma*x[t - 1] + eps_x[t]
        }
        # y = beta*[0x(1:T-1)] + eps_y
      }
    }
  }
}

# Compute statistics
# x_lag = cbind(x(1, x(1:T-1)
# x_lag_demeaned = x_lag - mean(x_lag)
# vhat_x = x - x(1)
# vhat_x_lag = [0vhat_x(1:T-1)]
#
# % Estimate Omega
# % True values: omegahat_yy = 1 omegahat_xx = 1 omegahat_xy = rho
#
# epshat_y = y - mean(y) - x_lag_demeaned.*(x_lag_demeaned'*y./(x_lag_demeaned'*x_lag_demeaned))
# epshat_x = vhat_x - vhat_x_lag.*(vhat_x_lag'*vhat_x./(vhat_x_lag'*vhat_x_lag))
#
# omegahat_yy = (epshat_y'*epshat_y)./T
#                omegahat_xy = (epshat_x'*epshat_y)./T
# omegahat_xx = (epshat_x'*epshat_x)./T
#
#                omegahat_yy_x = omegahat_yy - (omegahat_xy.^2)./omegahat_xx
#                rhohat = omegahat_xy ./ sqrt(omegahat_xx.*omegahat_yy)
#                rhobarhat = rhohat./sqrt(1 - rhohat.^2)
#
#                % Construct minimal sufficient statistics
#
#                Rhat_b = ((x_lag_demeaned'*(y - beta_0*x_lag)) ./ T - (omegahat_xy./2) .* ((vhat_x(T).^2)./(T.*omegahat_xx) - 1) ...
# + (omegahat_xy ./ ((T.^2).*omegahat_xx)) .* vhat_x(T) .* sum(vhat_x(1:T-1))) ./ sqrt(omegahat_xx.*omegahat_yy_x)
# Rhat_g = ((vhat_x(T).^2)./(T.*omegahat_xx) - 1) ./ 2 - rhobarhat .* Rhat_b
# Rhat_bb = (x_lag_demeaned'*x_lag_demeaned)./((T.^2).*omegahat_xx)
#            Rhat_gg = (vhat_x_lag'*vhat_x_lag)./((T.^2).*omegahat_xx)
#
# % Compute conditional p-values
#
# if rhobarhat < 0
#
# rtilde_b_threshold = 1 - exp(-(2.*Rhat_g+2.*rhobarhat.*Rhat_b+1).^(1./k))
#
# Zhat_bb = sqrt(Rhat_gg-Rhat_bb)
#
# Integral_tolerance = Integral_tolerance_initial
# Integral_iteration = 1
# Integral_full_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_full_plus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_full_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_reliability_plus(i,b_index) = (Integral_full_plus(i,b_index)>(1.e3.*Integral_tolerance))
#
# Integral_iteration = 1
# Integral_lower_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,rtilde_b_threshold,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_lower_plus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_lower_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,rtilde_b_threshold,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_tolerance = Integral_tolerance_initial
# Integral_iteration = 1
# Integral_full_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_full_minus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_full_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_reliability_minus(i,b_index) = (Integral_full_minus(i,b_index)>(1.e3.*Integral_tolerance))
#
# Integral_iteration = 1
# Integral_lower_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,rtilde_b_threshold,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_lower_minus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_lower_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,rtilde_b_threshold,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_lower(i,b_index) = Integral_lower_minus(i,b_index) + Integral_lower_plus(i,b_index)
# Integral_full(i,b_index) = Integral_full_minus(i,b_index) + Integral_full_plus(i,b_index)
#
# Integral_reliability(i,b_index) = Integral_reliability_plus(i,b_index).*Integral_reliability_minus(i,b_index)
#
# pvalue(i,b_index) = 1 - Integral_lower(i,b_index) ./ Integral_full(i,b_index)
#
# elseif rhobarhat == 0
#
# pvalue(i,b_index) = 1 - normcdf(Rhat_b./sqrt(Rhat_bb))
#
# elseif rhobarhat > 0
#
# rtilde_b_threshold = 1 - exp(-(2.*Rhat_g+2.*rhobarhat.*Rhat_b+1).^(1./k))
# rtilde_b_thresholds(i,b_index) = rtilde_b_threshold
#
# Zhat_bb = sqrt(Rhat_gg-Rhat_bb)
#
# Integral_tolerance = Integral_tolerance_initial
# Integral_iteration = 1
# Integral_full_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_full_plus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_full_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_reliability_plus(i,b_index) = (Integral_full_plus(i,b_index)>(1.e3.*Integral_tolerance))
#
# Integral_iteration = 1
# Integral_lower_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,Integral_lower_truncation_rtilde,rtilde_b_threshold,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_lower_plus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_lower_plus(i,b_index) = doublequad(@pvalue_integrand_plus_internal,Integral_lower_truncation_rtilde,rtilde_b_threshold,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_tolerance = Integral_tolerance_initial
# Integral_iteration = 1
# Integral_full_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_full_minus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_full_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,Integral_lower_truncation_rtilde,1-Integral_upper_truncation_rtilde,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_reliability_minus(i,b_index) = (Integral_full_minus(i,b_index)>(1.e3.*Integral_tolerance))
#
# Integral_iteration = 1
# Integral_lower_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,Integral_lower_truncation_rtilde,rtilde_b_threshold,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# while ((Integral_lower_minus(i,b_index)>=1.e3.*Integral_tolerance) + (Integral_iteration>=Adjustments_max)) == 0,
#
# Integral_iteration = Integral_iteration + 1
# Integral_tolerance = Integral_tolerance./Integral_tolerance_correction
# Integral_lower_minus(i,b_index) = doublequad(@pvalue_integrand_minus_internal,Integral_lower_truncation_rtilde,rtilde_b_threshold,Integral_lower_truncation_ttilde,1-Integral_upper_truncation_ttilde,Integral_tolerance,@quadl,Rhat_g,Rhat_bb,Rhat_gg,Zhat_bb,rhobarhat,k)
#
# end % while
#
# Integral_lower(i,b_index) = Integral_lower_minus(i,b_index) + Integral_lower_plus(i,b_index)
# Integral_full(i,b_index) = Integral_full_minus(i,b_index) + Integral_full_plus(i,b_index)
#
# Integral_reliability(i,b_index) = Integral_reliability_plus(i,b_index).*Integral_reliability_minus(i,b_index)
#
# pvalue[i, b_index] = 1 - Integral_lower[i,b_index] / Integral_full[i, b_index]
#
# end #% if rhobarhat
#
# end #% for i
#
# end #% for b_index
#
# Rejection_rate = rowSums(pvalue<=eta)/n
# Rejection_rate_trimmed = 1 - rowSums((pvalue.*Integral_reliability )>= eta)/sum(Integral_reliability)
# Trimming = 1 - rowSums(Integral_reliability)/n
#
# # % Compute equal tails two-sided p-values
#
# pvalue_equaltails_twosided <-  2*min(pvalue, 1-pvalue)
# Rejection_rate_equaltails_twosided  <-  rowSums(pvalue_equaltails_twosided <= eta)/n
#
# end #% for c_index
#
# end #% for rho_index

function f_plus = pvalue_integrand_plus_internal(rtilde_b,ttilde_gg,r_g,r_bb,r_gg,z_bb,rhobarhat,k)

tau = sqrt(-log(1-max(ttilde_gg,1.e-8)));
    
A(1,1) = (1 ./ tau) .* (sinh(2.*tau) + sin(2.*tau)) ./ (cosh(2.*tau) + cos(2.*tau));
A(1,2) = (1 ./ (tau.^2)) .* (2.*sinh(tau).*sin(tau)) ./ (cosh(2.*tau) + cos(2.*tau));
A(2,1) = A(1,2);
A(2,2) = (1 ./ (2.*(tau.^3))) .* (sinh(2.*tau) - sin(2.*tau)) ./ (cosh(2.*tau) + cos(2.*tau));
   
B(1,1) = (1 ./ tau) .* (sinh(2.*tau) - sin(2.*tau)) ./ (cosh(2.*tau) + cos(2.*tau));
B(1,2) = (1 ./ (tau.^2)) .* (1 - (2.*cosh(tau).*cos(tau)) ./ (cosh(2.*tau) + cos(2.*tau)));
B(2,1) = B(1,2);
B(2,2) = (1 ./ (tau.^2)) .* (1 - (1 ./ (2.*tau)) .* (sinh(2.*tau) + sin(2.*tau)) ./ (cosh(2.*tau) + cos(2.*tau)));

B_inv = inv(B);

ReMatrix = - inv(A*B_inv*A + B)*A*B_inv;
ImMatrix = B_inv + B_inv*A*ReMatrix;

ReMatrix_11 = ReMatrix(1,1);
ReMatrix_1221 = ReMatrix(1,2)+ReMatrix(2,1);
ReMatrix_22 = ReMatrix(2,2);

ImMatrix_11 = ImMatrix(1,1);
ImMatrix_1221 = ImMatrix(1,2)+ImMatrix(2,1);
ImMatrix_22 = ImMatrix(2,2);

Multiplier = 1 ./ (sign(cos(tau./2)) .* sqrt(cosh(sqrt(-2.*i.*(tau.^2)))) .* sqrt(det(A+i*B)));

Multiplier_real = real(Multiplier);
Multiplier_imag = imag(Multiplier);

r_b = ((-log(1-rtilde_b)).^k-2.*r_g-1)./(2.*rhobarhat);
q_g = (-log(1-rtilde_b)).^k;

exponent_plus_real = -(r_b.^2)./(2.*r_bb) + ReMatrix_11.*q_g + ReMatrix_1221.*sqrt(q_g).*z_bb + ReMatrix_22.*(z_bb.^2);

exponent_plus_imag = log(1-ttilde_gg).*r_gg + ImMatrix_11.*q_g + ImMatrix_1221.*sqrt(q_g).*z_bb + ImMatrix_22.*(z_bb.^2);

f_plus = exp(exponent_plus_real) .* (Multiplier_real .* cos(exponent_plus_imag) - Multiplier_imag .* sin(exponent_plus_imag)) ...
    .* ((-log(1-rtilde_b)).^(k./2-1) ./ ((1-rtilde_b).*(1-ttilde_gg))) ./ (abs(rhobarhat).*sqrt(r_bb).*z_bb);
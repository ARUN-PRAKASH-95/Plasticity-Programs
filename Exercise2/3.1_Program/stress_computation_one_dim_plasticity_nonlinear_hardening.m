function [sigma_n1,epsilon_pl_n1,alpha_n1,X_n1] = stress_computation_one_dim_plasticity_nonlinear_hardening(epsilon_n1,epsilon_pl_n,alpha_n,X_n,Emod,sig0,H,c,h,eta,DeltaY)
%
% computes the stress and the internal variables (plastic strain and
% isotropic hardening variable) at time t_n1 based on euler-implicit time
% integration
%
% compute trial stress
sigma_tr = Emod*(epsilon_n1-epsilon_pl_n);
X_n = H * epsilon_pl_n
R = h * alpha_n + DeltaY*(1-exp(-eta*alpha_n))
%
% increase in yield stress at t_n
[R_n,~] = nonlinear_isotropic_hardening_exp_type(alpha_n,DeltaY,eta,h,'funct');
%
% evaluate yield function with trial stress state
Pi_trial = abs(sigma_tr-X_n) - sig0 - R_n

% add your source code for the stress difference and the yield function at 
% the trial elastic state here
xi_trial = abs(sigma_tr-X_n) 


% check violation of yield surface at elastic trial state / corrector if
% necessary
if (Pi_trial<0)
  sigma_n1=sigma_tr;
  epsilon_pl_n1=epsilon_pl_n;
  alpha_n1=alpha_n;
  X_n1=X_n;
  
% add your source code for the elastic predictor / plastic corrector
% algorithm here, including a Newton-Raphson scheme for the determination
% of the incremental lagrange multiplier
else
  del_lambda = 0
  dif = 1
  while dif >= 10^-4
    
    del_lambda_f = xi_trial - sig0 - (h * alpha_n) - DeltaY*(1 - exp(-eta*alpha_n-eta*del_lambda)) / (Emod+H+h)
    del_lambda_der = eta*exp(-eta*del_lambda-eta*alpha_n)*DeltaY + (Emod+H+h)
    dif = del_lambda / del_lambda_der
    del_lambda_n1 = del_lambda_f - dif
   end
  epsilon_pl_n1 = epsilon_pl_n + del_lambda_n1*sign(sigma_tr - X_n)  % Implicit Euler for strain update
  alpha_n1 = alpha_n + del_lambda_n1                
  sigma_n1 = Emod*(epsilon_n1-epsilon_pl_n1) 
  X_n1 = H*(epsilon_pl_n1)
end


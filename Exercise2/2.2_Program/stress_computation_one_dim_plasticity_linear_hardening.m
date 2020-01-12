function [sigma_n1,epsilon_pl_n1,alpha_n1] = stress_computation_one_dim_plasticity_linear_hardening(epsilon_n1,epsilon_pl_n,alpha_n,Emod,sig0,H,h)
%
% computes the stress and the internal variables (plastic strain and
% isotropic hardening variable) at time t_n1 based on euler-implicit time
% integration
%
% compute trial stress
sigma_trial =  Emod*(epsilon_n1-epsilon_pl_n)
 % back-stress and increase in yield stress at t_n
X_n = H*epsilon_pl_n            % Back stress
R = h*alpha_n                   % Increase in elastic domain
% evaluate yield function with trial stress state
Pi_trial = abs(sigma_trial - X_n) - sig0 - R
% check violation of yield surface at elastic trial state / corrector if
% necessary
if (Pi_trial<0)
  sigma_n1= sigma_trial;
  epsilon_pl_n1=epsilon_pl_n;
  alpha_n1=alpha_n
else
  del_lambda = Pi_trial/(Emod+H+h)                     %plastic corrector
  epsilon_pl_n1 = epsilon_pl_n + del_lambda*sign(sigma_trial - X_n)  % Implicit Euler for strain update
  alpha_n1 = alpha_n + del_lambda                
  sigma_n1 = Emod*(epsilon_n1-epsilon_pl_n1)                         %Updated stress
end
% dummy code to provide functionality of the function (linear elastic material)

%



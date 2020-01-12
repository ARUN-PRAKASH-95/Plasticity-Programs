function [residual,jacobian] = residual_jacobian_yield_cond_nonlinear_hardening(delta_lambda,alpha_n,abs_xi_trial,abs_X_n,fact,Emod,sig0,H,c,h,eta,DeltaY)
%
% residual and corresponding jacobian for yield function \Phi_{n+1}=0
% determination of incremental lagrange multiplier
%
one_p_c_dlambda=1+c*delta_lambda;
alpha_n1_t=alpha_n+delta_lambda;
[R,dR_dalpha] = nonlinear_isotropic_hardening_exp_type(alpha_n1_t,DeltaY,eta,h,'both');
%
residual=abs_xi_trial+fact*(c*delta_lambda*abs_X_n)/one_p_c_dlambda-delta_lambda*(Emod+H/one_p_c_dlambda)-sig0-R;
jacobian=fact*(abs_X_n*c*one_p_c_dlambda-delta_lambda*abs_X_n*c^2)/one_p_c_dlambda^2 -Emod -(H*one_p_c_dlambda-H*delta_lambda*c)/one_p_c_dlambda^2- dR_dalpha;
%
end


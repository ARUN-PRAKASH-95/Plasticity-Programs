function [epsilon_total,sigma,X] = analyt_sol_ramp_strain_nonlinear_iso_lin_kin_hard_plast(n_mult,sig0,Emod,H,h,DeltaY,eta)
%
% compute the reference solution for a nonlinear isotropic and linear kinematic 
% hardening, one-dimensional plasticity model subject to monotone
% increasing strain history (ramp)
%
n_step_pl=10;
sign_sigma_m_x=1;
%
dim_less_time=1/n_mult:(1-1/n_mult)/(n_step_pl-1):1;
lambda_hist=zeros(1,n_step_pl);
residual_lambda=@(lambda,dl_time,sig0,Emod,n_mult,H,h) (Emod+H+h)*lambda+DeltaY*(1-exp(-eta*lambda))-sign_sigma_m_x*n_mult*sig0*(dl_time-1/n_mult);
options = optimset('Display','off');
for i=2:n_step_pl
%    residual_temp=residual_lambda(lambda_hist(i),dim_less_time(i),sig0,Emod,n,H,h,c);
    [lambda_temp,fval,exitflag]=fzero(@(lambda) residual_lambda(lambda,dim_less_time(i),sig0,Emod,n_mult,H,h),[0,n_mult*sig0/Emod],options);
    lambda_hist(i)=lambda_temp;
end
epsilon_plast=[0.0,lambda_hist*sign_sigma_m_x];
epsilon_total=[0.0,dim_less_time*n_mult*sig0/Emod];
%
sigma=Emod*(epsilon_total-epsilon_plast);
X=H*epsilon_plast;
end


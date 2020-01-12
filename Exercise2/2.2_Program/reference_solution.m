function [stress,epsilon_total] = reference_solution(T,time_prop,n_mult,Emod,sig0,H,h)
%
% computes the reference solution for strain history in terms of triangle
% wave function
%
t_star=1.0/n_mult;
t_d_star=1.0/n_mult*((2.0*n_mult+1.0+h*n_mult/(Emod+H+h)*(1.0-1.0/n_mult))-(n_mult*Emod)/(Emod+H+h)*(1.0-1.0/n_mult)*(1.0+H/Emod));
dim_less_time=[0.0,t_star,1.0,t_d_star,3.0];
epsilon_total=zeros(5,1);
for i=1:5
    epsilon_total(i)=comput_strain_prscr_comp(dim_less_time(i)*T,time_prop);
end
%
lagr_mp_T=n_mult*sig0/(Emod+H+h)*(1.0-1.0/n_mult);
lagr_mp_3T=lagr_mp_T+n_mult*sig0/(Emod+H+h)*(3.0-t_d_star);
lagr_multiplier=[0.0,0.0,lagr_mp_T,lagr_mp_T,lagr_mp_3T];
%
sign_sigma_backstress=[0.0,1.0,1.0,-1.0,-1.0];
epsilon_plastic_t=zeros(5,1);
for i=1:4
    epsilon_plastic_t(i+1)=epsilon_plastic_t(i)+sign_sigma_backstress(i+1)*(lagr_multiplier(i+1)-lagr_multiplier(i));
end
%
stress=Emod*(epsilon_total-epsilon_plastic_t);
%
end


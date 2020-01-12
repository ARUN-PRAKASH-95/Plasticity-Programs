clear all;
% *************************************************************************
% one-dimensional plasticity with nonlinear isotropic and kinematic hardening
% *************************************************************************
%
% material parameter
Emod=20000;           % Young's modulus in MPa
sig0=200;             % initial yield stress in MPa
H=12000;              % kinematic hardening modulus in MPa
h=10;                 % isotropic hardening modulus in MPa
eta=500;              % dimensionless shape parameter for nonlinear isotropic hardening
DeltaY=100;           % asymptotic increase in yield stress (nonlinear isotropic hardening)
c=0;                  % parameter for AF-kinematic hardening (not used)
%
% simulation parameters
T=10;                                 % time for load segment
n_cycle=1;                            % number of cycles
T_sim=4*n_cycle*T;                    % simulation time 
n_incr_p_cycle=120;                    % number of increments per cycle (15 is of interest)
n_incr_total=n_cycle*n_incr_p_cycle;  % total number of increments
n_mult=3;                             % multiple of the elastic strain at initial yield
strain_ampl=n_mult*sig0/Emod;         % strain amplitude
%
%% computations
%
dt=T_sim/n_incr_total;
t_n=0;
alpha=zeros(1,n_incr_total+1);
epsilon_pl=zeros(1,n_incr_total+1);
X_back_stress=zeros(1,n_incr_total+1);
epsilon_total=zeros(1,n_incr_total+1);
time=zeros(1,n_incr_total+1);
%
time_prop = loading_history(T,strain_ampl,n_cycle);
%
for i=1:n_incr_total
    t_n1=t_n+dt;
    time(i+1)=t_n1;
%
% assign variables from previous time step
    alpha_n=alpha(i);
    epsilon_pl_n=epsilon_pl(i);
    X_n=X_back_stress(i);
%
% get the strain value at t_n1
    epsilon_n1=comput_strain_prscr_comp(t_n1,time_prop);
    epsilon_total(i+1)=epsilon_n1;
%
% compute stress and internal variables
    [sigma_n1,epsilon_pl_n1,alpha_n1,X_n1] = ...
        stress_computation_one_dim_plasticity_nonlinear_hardening(epsilon_n1,...
        epsilon_pl_n,alpha_n,X_n,Emod,sig0,H,c,h,eta,DeltaY);
%
% store internal variables
    alpha(i+1)=alpha_n1;
    epsilon_pl(i+1)=epsilon_pl_n1;
    X_back_stress(i+1)=X_n1;
%
    t_n=t_n1;
%
end
%
sigma=Emod*(epsilon_total-epsilon_pl);
%
% compute reference solutions
% nonlinear isotropic, linear kinematic
[epsilon_total_ref_isotr_hard,sigma_ref_isotr_hard,X_ref_isotr_hard] = ...
    analyt_sol_ramp_strain_nonlinear_iso_lin_kin_hard_plast(n_mult,sig0,...
    Emod,H,h,DeltaY,eta);
%
figure(1)
clf;
plot(epsilon_total,sigma,'bo--')
hold on
plot(epsilon_total,X_back_stress,'r-')
xlabel('\epsilon')
ylabel('\sigma, X in MPa')
title('stress-strain curve')
%
figure(2)
clf;
plot(time/T,epsilon_total,'bx-')
grid on
xlabel('t/T')
ylabel('\epsilon')
title('strain history')
%
figure(3)
clf;
[R_tt,~] = nonlinear_isotropic_hardening_exp_type(alpha,DeltaY,eta,h,'funct');
upper_bound_elastic_domain=X_back_stress+sig0+R_tt;
lower_bound_elastic_domain=X_back_stress-sig0-R_tt;
ha = shadedplot(time/T, upper_bound_elastic_domain, lower_bound_elastic_domain, [0.85 0.85 0.85], 'b');
hold on
plot(time/T,X_back_stress,'r-')
xlabel('t/T')
ylabel('\sigma, X in MPa')
title('evolution of elastic domain')
%
figure(4)
clf;
alpha_t=0:4*strain_ampl/99:4*strain_ampl;
[R_t,~] = nonlinear_isotropic_hardening_exp_type(alpha_t,DeltaY,eta,h,'funct');
plot(alpha_t,sig0+R_t,'k-')
xlabel('\alpha')
ylabel('\sigma^{y}')
%
figure(5)
clf;
incr_oi=ceil(T/dt)+1;
if time(incr_oi)-T>sqrt(eps)
    incr_oi=incr_oi-1;
end
plot(epsilon_total(1:incr_oi),sigma(1:incr_oi),'b--')
hold on
plot(epsilon_total(1:incr_oi),X_back_stress(1:incr_oi),'r-')
plot(epsilon_total_ref_isotr_hard,sigma_ref_isotr_hard,'k*')
plot(epsilon_total_ref_isotr_hard,X_ref_isotr_hard,'g*')
xlabel('\epsilon')
ylabel('\sigma, X in MPa')
title('stress-strain curve')
%
legend('\sigma-\epsilon','X-\epsilon','\sigma-\epsilon (Ref., nl. iso, lin. kin)',...
    'X-\epsilon (Ref. nl. iso, lin. kin)','Location','northwest')
clear all;
% *************************************************************************
% one-dimensional plasticity with linear isotropic and kinematic hardening
% *************************************************************************
%
% material parameter
Emod=20000;           % Young's modulus in MPa
sig0=200;             % initial yield stress in MPa
H=500;                % kinematic hardening modulus in MPa
h=500;                % isotropic hardening modulus in MPa
%
% simulation parameters
T=10;                          % time for load segment
T_sim=4*T;                     % simulation time 
n_incr=150;                    % number of increments (15 is of interest)
n_mult=3;                      % multiple of the elastic strain at initial yield
strain_ampl=n_mult*sig0/Emod;  % strain amplitude
%
%% computations
%
dt=T_sim/n_incr;
t_n=0;
alpha=zeros(1,n_incr+1);
epsilon_pl=zeros(1,n_incr+1);
epsilon_total=zeros(1,n_incr+1);
time=zeros(1,n_incr+1);
%
time_prop = loading_history(T,strain_ampl);
%
for i=1:n_incr
    t_n1=t_n+dt;
    time(i+1)=t_n1;
%
% assign variables from previous time step
    alpha_n=alpha(i);
    epsilon_pl_n=epsilon_pl(i);
%
% get the strain value at t_n1
    epsilon_n1=comput_strain_prscr_comp(t_n1,time_prop);
    epsilon_total(i+1)=epsilon_n1;
%
% compute stress and internal variables
    [sigma_n1,epsilon_pl_n1,alpha_n1] = ...
        stress_computation_one_dim_plasticity_linear_hardening(epsilon_n1,...
        epsilon_pl_n,alpha_n,Emod,sig0,H,h);
%
% store internal variables
    alpha(i+1)=alpha_n1;
    epsilon_pl(i+1)=epsilon_pl_n1;
%
    t_n=t_n1;
%
end
%
sigma=Emod*(epsilon_total-epsilon_pl);
X_back_stress=H*epsilon_pl;
%
[stress_ref,epsilon_ref] = reference_solution(T,time_prop,n_mult,Emod,sig0,H,h);
%
figure(1)
clf;
plot(epsilon_total,sigma,'bo--')
hold on
plot(epsilon_ref,stress_ref,'k:')
plot(epsilon_total,X_back_stress,'r-')
xlabel('\epsilon')
ylabel('\sigma, X in MPa')
title('stress-strain curve')
%
figure(2)
clf;
plot(time/T,epsilon_total,'bx-')
xlabel('t/T')
ylabel('\epsilon')
title('strain history')
%
figure(3)
clf;
upper_bound_elastic_domain=X_back_stress+sig0+h*alpha;
lower_bound_elastic_domain=X_back_stress-sig0-h*alpha;
ha = shadedplot(time/T, upper_bound_elastic_domain, lower_bound_elastic_domain, [0.85 0.85 0.85], 'b');
hold on
plot(time/T,X_back_stress,'r-')
xlabel('t/T')
ylabel('\sigma, X in MPa')
title('evolution of elastic domain')
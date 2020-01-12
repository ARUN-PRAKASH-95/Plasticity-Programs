clear all;
% *************************************************************************
% test suite for three-dimensional von Mises plasticity 
% isochoric tension (fully-deformation controlled test, ramp function)
% *************************************************************************
%
% strain amplitude in terms of multiple of normalized yield stress
% strain_ampl=sigma_y0/2/mu*n_ampl
n_ampl=2;      % n_ampl>1
%
mat_param = inputmat();
xE = mat_param(1); xnu = mat_param(2); sigma_y0 = mat_param(3);
mu = xE/(2*(1+xnu));
%
% define loading
% 1: linear ramping of load
l11type=1;
l12type=3;
dt=0.01;
if l11type==1
    t11=[0 10];
    lam11=[0 n_ampl*sigma_y0/2/mu];
end
if l12type==3
    t12=[0 5 10];
    lam12=[0 0.0 0.01];
end
%
% computation of tangent moduli
ttype = 0; % 0: analytical
%
% path to auxiliary functions
addpath('tensor/');
%% computation
%
% dt=dt*10; t11=t11*10; t12=t12*10;
% prescribed load/time step
% start and end-time of loading, time-scale, no. of steps
ta=t11(1);
te=t11(end);
time=ta:dt:te;
steps=size(time,2)-1;
e11=loading(l11type,dt,t11,lam11);
e12=loading(l12type,dt,t12,lam12);
%
% initialize internal variables 
sdv=zeros(26,steps+1);
%
% initialise quantities for post-processing
s11=zeros(1,steps+1); s22=zeros(1,steps+1); s33=zeros(1,steps+1); s12=zeros(1,steps+1);
%
% initialize waitbar
wb=waitbar(0,'computation in progress...');

for n=1:steps
%
% display waitbar
    waitbar(n/steps);
% display current time step
    disp(['n = ', num2str(n)]);
%
    epsilon=[1;-0.5;-0.5;0;0;0]*e11(n+1)+[0;0;0;1;0;0]*e12(n+1);
%
% constitutive law: algorithmic stresses and moduli 
    [s,A,sdvup]=vmises_perzyna(epsilon,sdv(:,n),dt,ttype);
%
% update of internal variables after obtaining convergence
    sdv(:,n+1) = sdvup;
%    
% store quantities for post-processing
    s11(n+1)=s(1); s22(n+1)=s(2); s33(n+1)=s(3); s12(n+1)=s(4);
%
end
close(wb);
%
% open the file with write permission
fid = fopen('data-isochoric-tension-shear-test-perzyna.txt', 'w');
fprintf(fid, '%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n', [time;e11;e12;s11;s22;s33;s12]);
fclose(fid);
%% visualization
figure(1)
clf;
hold on
plot(e11,s11,'r-')
plot(e11,s22,'b-')
plot(e11,s33,'g-')
plot(e11,s12,'m-')
%
xlabel('\epsilon_{11}','FontSize',12)
ylabel('\sigma_{11}, \sigma_{22}, \sigma_{33} in MPa','FontSize',12)
legend('\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}','Location','East')
%
figure(2)
clf;
hold on
plot(e11,sdv(1,:),'r-')
plot(e11,sdv(2,:),'b-')
plot(e11,sdv(3,:),'g-')
plot(e11,sdv(4,:),'m-')
%
xlabel('\epsilon_{11}','FontSize',12)
ylabel('\epsilon_{11}^{pl}, \epsilon_{22}^{pl}, \epsilon_{33}^{pl}, \epsilon_{12}^{pl}','FontSize',12)
legend('\epsilon_{11}^{pl}','\epsilon_{22}^{pl}','\epsilon_{33}^{pl}','\epsilon_{12}^{pl}','Location','NorthWest')
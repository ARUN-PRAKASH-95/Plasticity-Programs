function [sig6,A66,sdvl]=vmises_perzyna(eps6,sdvl,dt,ttype)

%==========================================================================
% vmises.m (version 4)
%
% standard radial return algorithm for small strain plasticity with
% linear isotropic and kinematic hardening
%==========================================================================
% 
% coded by: B. Kiefer 19 Feb 2013
%           last modification: 19 Feb 2013
%
% comments: -> execution:
%              drive
% 
%==========================================================================
tol=1e-8;

% material parameters
matp = inputmat();
xE        = matp(1);
xnu       = matp(2);
xsigy0    = matp(3);
xH        = matp(4);
xh        = matp(5);
eta       = matp(6);
xmu  = xE/(2*(1+xnu));
xk   = xE/(3*(1-2*xnu));
%
% characteristic time
tau=eta/2/xmu;
%
% general 
ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
xid = eye(3);

% restore the strain tensor
eps = [eps6(1) eps6(4)/2 eps6(6)/2;
       eps6(4)/2 eps6(2) eps6(5)/2;
       eps6(6)/2 eps6(5)/2 eps6(3)];

% restore the internal variables at tn
epspn   = [sdvl(1) sdvl(4)/2 sdvl(6)/2;
           sdvl(4)/2 sdvl(2) sdvl(5)/2;
           sdvl(6)/2 sdvl(5)/2 sdvl(3)];
Balphan = [sdvl(1+6) sdvl(4+6)/2 sdvl(6+6)/2;
           sdvl(4+6)/2 sdvl(2+6) sdvl(5+6)/2;
           sdvl(6+6)/2 sdvl(5+6)/2 sdvl(3+6)];
alphan  =  sdvl(13);
%
epsp = epspn;              % plastic strain tensor
Balpha = Balphan;          % tensorial internal variable alpha
alpha = alphan;            % scalar hardening variable alpha
%
% linear elastic behavior
deps = eps - 1/3*trace(eps)*xid;
%
% % TRIAL STEP
% deps = ................ deviatoric part of the strain tensor
deps = eps - 1/3*trace(eps)*xid;
% dsigtr = .............. deviatoric part of the trial stress tensor
dsigtr = 2*xmu*(deps - epsp);
% Bbetatr = ..............trial value of the back-stress tensor
Bbetatr = xH * Balpha;
% betatr = ...............trial value of the increase in yield stress (scalar beta)
betatr = xh * alpha;
% xitr = .................trial value of deviatoric stress difference
xitr = dsigtr - Bbetatr;
% nxitr = ................norm of xitr (use function t2_contr_t2(A,B) contained in ./tensor/ to compute);
nxitr = sqrt(t2_contr_t2(xitr,xitr));
% ntr = ..................flow direction from trial state
ntr = xitr/nxitr;
% %
% % trial yield function  
% phitr = ................trial value of the yield function
phitr = nxitr - sqrt(2/3)*(xsigy0+betatr);
%
% decide if elastic or elastic-plastic step
if phitr<0
  dsig = dsigtr;
  C_dev = xk*t2_otimes_t2(xid,xid) + 2*xmu*getP4sym();
  epsp = epspn;              % plastic strain tensor
  Balpha = Balphan;          % tensorial internal variable alpha
  alpha = alphan;            % scalar hardening variable alpha
% restore stress tensor as vector
else
  l = tau;
  gamman = phitr/(2*xmu*l + 2*xmu + xH + (2/3)*xh);
  dsig = dsigtr - 2*xmu*gamman*ntr;
  epsp = epspn + gamman*ntr;
  Balpha = Balphan + gamman*ntr;
  alpha = alphan + sqrt(2/3)*gamman;
  a = phitr/nxitr;
  b = 1/(l+1+(xH/2*xmu)+(xh/3*xmu));
  beta1 = 1 - (a*b);
  beta2 = (1-a)*b;
  %
% for tangent computation use the auxiliary functions (contained in ./tensor)
%     t2_otimes_t2(A,B)
%     getP4sym()
%
  C_dev = 2*xmu*beta1*getP4sym() - 2*xmu*beta2*t2_otimes_t2(ntr,ntr);
endif
sig = dsig + xk*trace(eps)*xid;
C = C_dev + xk*t2_otimes_t2(xid,xid);


sig6=zeros(6,1);
for i=1:6
    sig6(i) = sig(ii(i),jj(i));
end

% restore stiffness tensor as matrix
ii = [1,2,3,1,2,1];
jj = [1,2,3,2,3,3];
A66=zeros(6,6);
%if ttype==0
    for i=1:6
        for j=1:6
        A66(i,j) = C(ii(i),jj(i),ii(j),jj(j));
        end
    end
%end

% store history variables
sdvl(1:6,1)=[epsp(1,1) epsp(2,2) epsp(3,3)...
             2*epsp(1,2) 2*epsp(2,3) 2*epsp(1,3)]';
sdvl(7:12,1)=[Balpha(1,1) Balpha(2,2) Balpha(3,3)...
             2*Balpha(1,2) 2*Balpha(2,3) 2*Balpha(1,3)]';
sdvl(13) = alpha;
end

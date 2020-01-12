function top = numtoploc2(eps6,sdvl,ttype,s)

%-------------------------------------------------------------------------
% returns tangent operator associated with Cauchy-stress
% (numerical (consistent) tangent operator conjugate to Cauchy-stress
% stress, for details: see Miehe, CMAME, 134, 223-240 (1996))
% for small strains
%-------------------------------------------------------------------------

%numt = 1;

% pertubation parameter
eps=1e-8;

% preallocation of matrices
top=zeros(6,6);

for i=1:6

    % 1. calculation of pertubed strain tensor (ep)
    if i==1
        eps6p=eps6+[eps; 0; 0; 0; 0; 0];
    elseif i==2
        eps6p=eps6+[0; eps; 0; 0; 0; 0];
    elseif i==3
        eps6p=eps6+[0; 0; eps; 0; 0; 0];
    elseif i==4
        eps6p=eps6+[0; 0; 0; 0.5*eps; 0; 0];
    elseif i==5
        eps6p=eps6+[0; 0; 0; 0; 0.5*eps; 0];
    elseif i==6
        eps6p=eps6+[0; 0; 0; 0; 0; 0.5*eps];
    end

    % 2. calculation of pertubed stresses (sigp, in Voigt-notation)
    %[pp,dummy1]=neoHooke(F9p,matp);
    [sp,dummy1,dummy2]=vmises(eps6p,sdvl,ttype);
    %[siggenp, dummy1, dummy2] = damage_nonlocal2(F9p, matparam,stateVar,numt);

    % 3. calculation of moduli
    testtop(:,i)=(sp(:,1)-s(:,1));
    top(:,i)=(sp(:,1)-s(:,1))/eps;

end

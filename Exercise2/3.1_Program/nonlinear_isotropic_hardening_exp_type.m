function [R,dR_dalpha] = nonlinear_isotropic_hardening_exp_type(alpha,DeltaY,eta,h,output)
%
% computes the increase of yield stress based on exponential type hardening
% law, with asymptotically linear hardening
%
switch output
    case 'funct'
        R=0; % add your code for a proper definition of the hardening function
        dR_dalpha=0;
    case 'deriv'
        R=0;
        dR_dalpha=0; % add your code for a proper definition of the derivative of the hardening function
    case 'both'
        % add your code for a proper definition of the hardening function
        % and its derivative
        R=0;
        dR_dalpha=0;
    otherwise
        error('unknown output identifier')
end
%
end


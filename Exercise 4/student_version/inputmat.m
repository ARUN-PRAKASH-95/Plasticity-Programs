function [matp] = inputmat()
%inputmat.m define material parameters

matp(1)  = 210.0e3;         % xE
matp(2)  = 0.33;            % xnu 
matp(3)  = 200.0;           % xsigy0
matp(4)  = 50.0*matp(3);    % xH
matp(5)  = 10.0*matp(3);    % xh

end

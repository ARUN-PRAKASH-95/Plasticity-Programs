function time_prop = loading_history(T,ampl)
%
% load history in terms of applied strain
%
% time_prop.........matrix containing in each row (time, prop-factor) value
%
time_prop=[[0.0,0.0];
           [T,ampl];
           [3*T,-ampl];
           [4*T,0.0]];
%
end


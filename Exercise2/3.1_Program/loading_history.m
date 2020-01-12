function time_prop = loading_history(T,ampl,n_cycle)
%
% load history in terms of applied strain
%
% time_prop.........matrix containing in each row (time, prop-factor) value
%
time_prop=[[0.0,0.0];
           [T,ampl];
           [3*T,-ampl];
           [4*T,0.0]]; % end of cycle one
for i=2:n_cycle
    time_prop(4+3*(i-2)+1,:)=[(4*(i-1)+1)*T,ampl];
    time_prop(4+3*(i-2)+2,:)=[(4*(i-1)+3)*T,-ampl];
    time_prop(4+3*(i-1),:)=[4*i*T,0.0];
end
%            [5*T,ampl];
%            [7*T,-ampl];
%            [8*T,0.0]; % end of cycle two
%            [9*T,ampl];
%            [11*T,-ampl];
%            [12*T,0.0]; % end of cycle three
%            [13*T,ampl];
%            [15*T,-ampl];
%            [16*T,0.0]]; % end of cycle four
%
end


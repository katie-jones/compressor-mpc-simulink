% Parameters of large exit tank 

function [Pout,Vtank,D,m_out_c] = const_tank()
%#eml

Pout = 1; % pressure at exit of tank

Vtank = 15 * pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 5.940; % volume, tank

% coefficients of the outlet valve characteristic
D = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];
% D = [1.5*[1 1 1 1], 1 1 1 1].*[-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];

m_out_c = 0.017;

end
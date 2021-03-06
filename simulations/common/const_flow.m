function [a,Pin,Ptank,V1,V2,AdivL] = const_flow()
%#eml

a = 340; % speed of sound

Pin = 1; % pressure at inlet
Ptank = -1; % pressure in tank -- shouldn't use

V1 = 2*pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 8.191; % volume of the inlet pipe, 6m
V2 =  pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 5.940; % volume, tank

AdivL =     pi * (0.08 / 2)^2 / 3 * 0.1;  % Duct cross section divided by length



end
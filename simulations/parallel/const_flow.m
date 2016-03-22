function [a,Pin,Ptank,V1,V2,AdivL,uoff1,uoff2] = const_flow()

a = 340; % speed of sound

Pin = 1; % pressure at inlet
Ptank = -1; % pressure in tank -- shouldn't use

V1 = pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 8.191; % volume of the inlet pipe, 6m
V2 =  pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 5.940; % volume, tank

AdivL =     pi * (0.08 / 2)^2 / 3 * 0.1;  % Duct cross section divided by length

ud1 = 0.393; % setting of discharge valve 1
ud2 = 0.393; % setting of discharge valve 2

uoff1 = [0.304, 0.405, ud1, 0, 0]'; % offset applied to calculated inputs
uoff2 = [0.304, 0.405, ud2, 0, 0]';

end
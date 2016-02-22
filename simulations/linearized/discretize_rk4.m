% Discretize system given by A,B,C using 4th order runge kutta
function [Ad,Bd,Cd] = discretize_rk4(A,B,C)
%#eml
Ts = 0.05;
xsize = 5;

A2 = A^2;
A3 = A^3;
A4 = A^4;

Ad = eye(xsize) + Ts*A + Ts^2/2*A2 + Ts^3/6*A3 + Ts^4/24*A4;
Bd = Ts*B + Ts^2/2*A*B + Ts^3/6*A2*B + Ts^4/24*A3*B;
Cd = C;

end
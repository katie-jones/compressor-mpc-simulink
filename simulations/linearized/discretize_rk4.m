% Discretize system given by A,B,C using 4th order runge kutta
function [Ad,Bd,Cd,fd] = discretize_rk4(A,B,C,f,Ts)
%#eml
% Ts = 0.05;
xsize = 5;

A2 = A*A;
A3 = A2*A;
A4 = A3*A;

Acom = Ts*eye(xsize) + Ts^2/2*A + Ts^3/6*A2 + Ts^4/24*A3;

% Ad = eye(xsize) + Ts*A + Ts^2/2*A2 + Ts^3/6*A3 + Ts^4/24*A4;
% Bd = Ts*B + Ts^2/2*A*B + Ts^3/6*A2*B + Ts^4/24*A3*B;
% Cd = C;
Ad = eye(xsize) + Acom*A;
Bd = Acom*B;
Cd = C;
fd = Acom*f(:);

end
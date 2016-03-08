% Discretize system given by A,B,C using 4th order runge kutta
function [Ad,Bd,Cd,fd] = discretize_rk4(A,B,C,f,Ts)
xsize = 5;

A2 = A*A;
A3 = A2*A;

Acom = Ts*eye(xsize) + Ts^2/2*A + Ts^3/6*A2 + Ts^4/24*A3;

Ad = eye(xsize) + Acom*A;
Bd = Acom*B;
Cd = C;
fd = Acom*f(:);

end
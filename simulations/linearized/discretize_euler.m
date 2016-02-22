% Discretize system given by A,B,C using 4th order runge kutta
function [Ad,Bd,Cd] = discretize_euler(A,B,C)
Ts = 0.05;
xsize = 5;

Ad = eye(xsize) + Ts*A;
Bd = Ts*B;
Cd = C;

end
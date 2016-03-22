function [Ts, n_delay,dsize,ucontrolsize,p,m,UWT,YWT] = const_mpc2()
Ts = 0.05;

n_delay = [0, 20];

dsize = 2;

ucontrolsize = 2;

p = 100;
m = 2;

UW = diag([100 1e4 100 1e4]');
YW = diag([0.1 1 0.1 1]');

UWT = kron(eye(m),UW);
YWT = kron(eye(p),YW);

end
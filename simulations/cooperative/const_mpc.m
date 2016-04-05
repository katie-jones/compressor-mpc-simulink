function [n_delay,dsize,ucontrolsize,p,m,UWT,YWT] = const_mpc()
%#eml

dsize = 2;

ucontrolsize = 2;

n_delay = [0, 40];

p = 50;
m = 1;

UW = [5e2 1e6];
YW = [1 1 0 0];

UWT = kron(eye(m),diag([UW]'));
YWT = kron(eye(p),diag(YW'));

end
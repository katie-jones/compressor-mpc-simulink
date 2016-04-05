function [n_delay,dsize,ucontrolsize,p,m,UWT,YWT] = const_mpc()
%#eml

dsize = 2;

ucontrolsize = 2;

n_delay = [0, 40];

p = 100;
m = 2;

UW = [2e3 8e4];
YW = [1 1 0.1 0.5];

UWT = kron(eye(m),diag([UW,UW]'));
YWT = kron(eye(p),diag(YW'));

end
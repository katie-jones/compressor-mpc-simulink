function [n_delay,dsize,ucontrolsize,p,m,UWT,YWT] = const_mpc()
%#eml

dsize = 2;

ucontrolsize = 2;

n_delay = [0, 40];

p = 100;
m = 2;

UW = [1e4 7e5];
YW = [1 1e-5 1e4];

UWT = kron(eye(m),diag([UW]'));
YWT = kron(eye(p),diag(YW'));

end
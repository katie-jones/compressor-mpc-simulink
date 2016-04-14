function [n_delay,dsize,ucontrolsize,p,m,UWT,YWT] = const_mpc()
%#eml

dsize = 2;

ucontrolsize = 2;

n_delay = [0, 40];

p = 100;
m = 2;

UW = [1e4 2e5];
YW = [1 1 1e-5 1e4];

UWT = kron(eye(m),diag(UW));
YWT = kron(eye(p),diag(YW'));
% UWT = kron([eye(m-1), zeros(m-1,1); zeros(1,m-1), 10],diag(UW'));
% YWT = kron(diag(1/10*logspace(1,0,p)),diag(YW));

end

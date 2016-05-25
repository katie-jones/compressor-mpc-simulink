% Constants for compressor simulation
function [orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m,UWT,YWT] = const_mpc()
%#eml
orig_xsize = 5;
ysize = 2;
dsize = 2;
usize = 2;

n_delay = [0; 5];
xsize = orig_xsize+dsize+n_delay(1)+n_delay(2);

Ts = 0.05;

p = 100;
m = 2;

% YW=diag([1 200]');
% UW=diag([200 400]');
% YW = diag([1e4 1e3]');
UW = diag([100 1e4]');
YW = diag([0.1 1]');
% UW = diag([1e3 1e5]');
% UW = diag([20 30]');
% YW = diag([100 1000]');


YWT = kron(eye(p),YW);
UWT = kron(eye(m),UW);



end

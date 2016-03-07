% Constants for compressor simulation
function [orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m,UWT,YWT] = constants()
%#eml
orig_xsize = 5;
ysize = 2;
dsize = 2;
usize = 2;

n_delay = [0; 20];
xsize = orig_xsize+dsize+n_delay(1)+n_delay(2);

Ts = 0.05;

p = 100;
m = 2;

% YW=diag([1e1 1e3]');
% UW=diag([50 1]');
% YW = diag([1e4 1e3]');
UW = diag([50 100]');
YW = diag([1e4 1e5]');
% UW = diag([1e3 1e5]');
% UW = diag([1 1]');
% YW = diag([1 1]');

YWT = kron(eye(p),YW);
UWT = kron(eye(m),UW);

% add weights for keeping last input for p-m iterations
% UWT(end-usize+1:end,end-usize+1:end) = (p-m+1)*UW;


end

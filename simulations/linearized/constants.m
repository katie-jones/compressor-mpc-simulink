% Constants for compressor simulation
function [orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m] = constants()
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
end
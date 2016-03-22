function [Ts,comp_xsize,xsize,total_xsize,comp_usize,comp_ysize] = mpc_constants()
Ts = 0.05;

comp_xsize = 5;
xsize = 2*comp_xsize + 1;
total_xsize = xsize;

comp_usize = 5;
comp_ysize = 2;

end
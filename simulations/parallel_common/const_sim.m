function [Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2, ud] = const_sim()
%#eml

Ts = 0.05;

xsize_comp = 5;
xsize = 2*xsize_comp + 1;

usize_comp = 5;

ysize_comp = 2;
ysize = 2*ysize_comp;

uoff1 = [0.304, 0.43, 1, 0, 0]'; % offset applied to calculated inputs
uoff2 = uoff1;

ud = 0.7;

end


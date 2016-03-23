function [xsize_comp, xsize, usize_comp, usize, ysize_comp, ysize, uoff1, uoff2, ud] = const_sim()
%#eml

xsize_comp = 5;
xsize = 2*xsize_comp + 1;

usize_comp = 5;
usize = 2*usize_comp + 1;

ysize_comp = 2;
ysize = 2*ysize_comp;

uoff1 = [0.304, 0.405, 0.55, 0, 0]'; % offset applied to calculated inputs
uoff2 = [0.304, 0.405, 0.55, 0, 0]';

ud = 0.9;

end


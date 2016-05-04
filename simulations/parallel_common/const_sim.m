function [Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2, ud] = const_sim()
%#eml

Ts = 0.05;

xsize_comp = 5;
xsize = 2*xsize_comp + 1;

usize_comp = 4;

ysize_comp = 2;
ysize = 2*ysize_comp;

uoff1 = [0.304, 0.405, 1, 0]'; % offset applied to calculated inputs
uoff2 = [0.304, 1, 0.393, 0]';

ud = uoff2(2);

end


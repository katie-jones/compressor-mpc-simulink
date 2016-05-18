% Setup for parallel compressor simulation
clear all
addpath('../common')

Td1 = [0 0;
    5 0];
Td2 = [0 0;
    5 0];
ur1 = [0 0];
ur2 = [0 0];

% linearization point
 
% xinit = [0.898; 1.12; 0.152; 438; 0; 
%     0.903; 1.128; 0.152; 438; 0;
%     P_D];
xinit = [  0.86733   1.0308  0.17648   394.96        0  0.99889   1.1872  0.17648   394.96        0 ]';


yref = [ 1.0425 7.6617 1.1784 7.6617 1.0116 ];

[Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2] = const_sim();


n_delay = [0; 40];
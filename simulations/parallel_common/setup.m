% Setup for parallel compressor simulation
clear all
addpath('../common')

Td1 = [0 0];
Td2 = [0 0];
ur1 = [0 0];
ur2 = [0 0];

% linearization point
P_D = 0.9925;
  
xinit = [0.898; 1.12; 0.152; 438; 0; 
    0.903; 1.128; 0.152; 438; 0;
    P_D];

[Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2] = const_sim();


n_delay = [0; 40];
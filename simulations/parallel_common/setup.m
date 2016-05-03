% Setup for parallel compressor simulation
clear all
addpath('../common')

Td1 = [0 0;
    5 0];
Td2 = [0 0];
ur1 = [0 0];
ur2 = [0 0];

% linearization point
  
xinit = [0.890; 1.096; 0.159; 426; 0; 
    0.928; 1.142; 0.159; 426; 0];

[Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2] = const_sim();

ud = 0.7;

P_D = 0.954;

n_delay = [0; 40];
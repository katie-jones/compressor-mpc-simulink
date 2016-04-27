% Setup for parallel compressor simulation
clear all
addpath('../common')

Td1 = [0 0];
Td2 = [0 0];
ur1 = [0 0];
ur2 = [0 0];

% linearization point
  
xinit = [0.887; 1.086; 0.161; 421; 0; 
    0.938; 1.148; 0.161; 421; 0];

[Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2] = const_sim();


ud = 0.7;

P_D = 0.954;

n_delay = [0; 40];
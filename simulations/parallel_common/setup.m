% Setup for parallel compressor simulation
clear all
addpath('../common')

Td = 0;
u_rec = 0.0; % recycle opening

P_D = 1.08; % initialize tank pressure

% linearization point
  
x_init_lin = [0.899; 1.126; 0.15; 440; 0];

[Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2, ud] = const_sim();

u_d = ud;

n_delay = [0; 40];
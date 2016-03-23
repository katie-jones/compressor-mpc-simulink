% Setup for parallel compressor simulation
clear all
global P_D 
addpath('../common')

Td = 0.0; % torque input
u_rec = 0.0; % recycle opening
u_d = 0.9; % discharge valve opening

P_D = 1.08; % initialize tank pressure

Ts = 0.05;

% linearization point
  
x_init_lin = [0.899; 1.126; 0.15; 440; 0];

u_out = 0.55;
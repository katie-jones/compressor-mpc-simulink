% Setup for parallel compressor simulation
clear all
global P_D 
addpath('../common')

Td = 0.1; % torque input
u_rec = 0.1; % recycle opening
u_d = 0.9; % discharge valve opening

P_D = 1.07; % initialize tank pressure

Ts = 0.05;

% linearization point
  
x_init_lin = [0.912; 1.17; 0.14; 465; 0];

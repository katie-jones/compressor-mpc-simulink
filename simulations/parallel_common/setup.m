% Setup for parallel compressor simulation
clear all
global P_D 
addpath('../common')

% Td = [zeros(100,1); 0.1]; % torque input
Td = 0;
u_rec = 0.0; % recycle opening
u_d = 0.7; % discharge valve opening

P_D = 1.08; % initialize tank pressure

Ts = 0.05;

% linearization point
  
x_init_lin = [0.899; 1.126; 0.15; 440; 0];

uoff1 = [0.304 0.43 1 0 0];
uoff2 = uoff1;


function [J,tau_r,A,C,m_in_c,m_rec_ss_c,D,m_out_c,T_ss_c,SD_c,torque_drive_c] = comp_coeffs()
%#eml

J          = (0.4 + 0.2070) *0.4; % 0.8;  % shaft inertia: compressor and rotor
tau_r =   1/0.5 * 1;     % time constant of the recycle valve

%%% coefficients of the compressor map 
A = [0.000299749505193654     -0.000171254191089237      3.57321648097597e-05    -9.1783572200945e-07 ...
    -0.252701086129365         0.136885752773673         -0.02642368327081        0.00161012740365743 ...
    54.8046725371143          -29.9550791497765          5.27827499839098         0.693826282579158];

% coefficients of the inlet valve characteristic
C = [-0.423884232813775         0.626400271518973       -0.0995040168384753      0.0201535563630318 ...
 -0.490814924104294         0.843580880467905        -0.423103455111209      0.0386841406482887];

% coefficients of the outlet valve characteristic
D = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];

m_in_c = 0.0051;

m_rec_ss_c = [0.0047, 0.0263];

m_out_c = 0.017;

T_ss_c = [2.5543945754982,  47.4222669576423, 0.6218];

SD_c = [5.55, 0.66];

torque_drive_c = 15000;
 
end
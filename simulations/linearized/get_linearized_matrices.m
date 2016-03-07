% Function to linearize compressor system about point xbar, ubar
% Have: 
%       xdot = f(x, u)
%       y = g(x)
% 
% Output: A, B, C such that
%       xdot = f(xbar, ubar) + A*(x-xbar) + B*(u-ubar)
%       y = C*(x-xbar)
%
% Output system is still continuous time
%
function [A, B, C] = get_linearized_matrices(x,u)
%#eml
Ts = 0.05;
p1 = x(1);
p2 = x(2);
m_comp = x(3);
omega_comp = x(4);
m_rec = x(5);

torque_drive_in = u(1)+0.304;
Recycle_opening = u(2);
Inflow_opening = 0.405;
Outflow_opening = 0.393;
% Inflow_opening = u(2);
% Outflow_opening = u(3); 
% Recycle_opening = u(4);

% torque_drive = torque_drive * 15000 / (2 * pi * 50);
% torque_drive = torque_drive * (2 * pi * 50 / omega_comp);

%-------------------------------------------------------------------------%
%----------------------------- CONSTANTS ---------------------------------%
%-------------------------------------------------------------------------%
SpeedSound = 340;     % speed of sound
In_pres    = 1.0;  % input pressure
Out_pres   = 1.0; % output pressure

VolumeT1 = 2*pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 8.191; % volume of the inlet pipe, 6m
VolumeT2 =  pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 5.940; % volume, tank
AdivL =     pi * (0.08 / 2)^2 / 3 * 0.1;  % Duct cross section divided by length

J          = (0.4 + 0.2070) *0.4; % 0.8;  % shaft inertia: compressor and rotor
tauRecycle =   1/0.5 * 1;     % time constant of the recycle valve

% coefficients of the compressor map 
A_coeff = [0.000299749505193654     -0.000171254191089237      3.57321648097597e-05    -9.1783572200945e-07 ...
    -0.252701086129365         0.136885752773673         -0.02642368327081        0.00161012740365743 ...
    54.8046725371143          -29.9550791497765          5.27827499839098         0.693826282579158];

% inflow valve coefficients
C_coeff = [-0.423884232813775         0.626400271518973       -0.0995040168384753      0.0201535563630318 ...
     -0.490814924104294         0.843580880467905        -0.423103455111209      0.0386841406482887];

% outflow valve coefficients
D2 = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];

%-------------------------------------------------------------------------%
%---------------------------- DERIVATIVES --------------------------------%
%-------------------------------------------------------------------------%

% Partial derivatives of p1
dp1dot = [ SpeedSound * SpeedSound / VolumeT1 * 1e-5 * (-1/2*sign(In_pres-p1)*100/sqrt(abs(In_pres*100-p1*100))) * (C_coeff(1)*Inflow_opening^3 + C_coeff(2)*Inflow_opening^2 + C_coeff(3)*Inflow_opening + C_coeff(4)),...
         0,...
         -SpeedSound * SpeedSound / VolumeT1 * 1e-5,...
         0,...
         SpeedSound * SpeedSound / VolumeT1 * 1e-5
         ];

% Partial derivatives of p2
dp2dot = [ 0,...
         -SpeedSound * SpeedSound / VolumeT2 * 1e-5 * (1/2*100/sqrt(abs(p2*100-Out_pres*100))) * (D2(1)*Outflow_opening^3 + D2(2)*Outflow_opening^2 + D2(3)*Outflow_opening + D2(4)),...
         SpeedSound * SpeedSound / VolumeT2 * 1e-5,...
         0,...
         -SpeedSound * SpeedSound / VolumeT2 * 1e-5
         ];
     

% Intermediate calculations
dM_dmcomp = [3*omega_comp^2*m_comp^2, 2*omega_comp^2*m_comp, omega_comp^2, 0,...
            3*omega_comp*m_comp^2, 2*omega_comp*m_comp, omega_comp, 0,...
            3*m_comp^2, 2*m_comp, 1, 0]';
dM_dwcomp = [2*omega_comp*m_comp^3, 2*omega_comp*m_comp^2, 2*omega_comp*m_comp, 2*omega_comp,...
            m_comp^3, m_comp^2, m_comp, 1,...
            0, 0, 0, 0]';
        
M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
    omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
    m_comp^3 m_comp^2 m_comp 1]';

p_ratio = A_coeff * M;


% Partial derivatives of qc
dqcdot = [ AdivL*(p_ratio*1e5),...
         -AdivL*1e5,...
         AdivL * (p1*1e5)*A_coeff*dM_dmcomp,...
         AdivL * (p1*1e5)*A_coeff*dM_dwcomp,...
         0
         ];

% Partial derivatives of omegac
domegacdot = [ 0, 0, -1.0/J * 47.4222669576423, -1.0/J*torque_drive_in*15000/omega_comp^2, 0];

% Partial derivatives of qr
dqrdot = [ -tauRecycle * (0.0047 * 1/2 * Recycle_opening / sqrt(p2*1e5 - p1*1e5) * 1e5),...
         tauRecycle * (0.0047 * 1/2 * Recycle_opening / sqrt(p2*1e5 - p1*1e5) * 1e5),...
         0, 0,...
         -tauRecycle
         ];
     
     
%-------------------------------------------------------------------------%
%------------------------------ OUTPUTS ----------------------------------%
%-------------------------------------------------------------------------%
% Linearized A matrix
A = [dp1dot;
    dp2dot;
    dqcdot;
    domegacdot;
    dqrdot];

% Linearized B matrix
B = [0, 0;
    0, 0;
    0, 0;
    1.0/J*15000/omega_comp, 0;
    0, tauRecycle*0.0047 * sqrt(p2*1e5 - p1*1e5)];

% Linearized C matrix
C = [0, 1, 0, 0, 0;
    p2/(5.55*p1^2), -1/(5.55*p1), 1, 0, 0];

C(2,:) = 100*C(2,:); % SD output is multiplied by factor 100



end
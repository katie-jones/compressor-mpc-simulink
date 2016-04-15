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

p1 = x(1);
p2 = x(2);
m_comp = x(3);
omega_comp = x(4);
% m_rec = x(5);

torque_drive_in = u(1);
% Inflow_opening = 0.405;
% Outflow_opening = 0.393;
Inflow_opening = u(2);
Outflow_opening = u(3); 
Recycle_opening = u(4);
dummy = u(5);

% torque_drive = torque_drive * 15000 / (2 * pi * 50);
% torque_drive = torque_drive * (2 * pi * 50 / omega_comp);

%-------------------------------------------------------------------------%
%----------------------------- CONSTANTS ---------------------------------%
%-------------------------------------------------------------------------%
[J,tauRecycle,A_coeff,C_coeff,~,m_rec_ss_c,D2,~,T_ss_c,SD_c,torque_drive_c] = comp_coeffs();
[SpeedSound,In_pres,Out_pres,VolumeT1,VolumeT2,AdivL] = const_flow();

if dummy > 0
    Out_pres = dummy;
end

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
domegacdot = [ 0, 0, -1.0/J * T_ss_c(2), -1.0/J*torque_drive_in*torque_drive_c/omega_comp^2, 0];

% Partial derivatives of qr
dqrdot = [ -tauRecycle * (m_rec_ss_c(1) * 1/2 * Recycle_opening / sqrt(p2*1e5 - p1*1e5) * 1e5),...
         tauRecycle * (m_rec_ss_c(1) * 1/2 * Recycle_opening / sqrt(p2*1e5 - p1*1e5) * 1e5),...
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
x0 = 1e-2; % deadzone cutoff
dmrur = tauRecycle * m_rec_ss_c(1) * sqrt(p2*1e5 - p1*1e5);
if Recycle_opening < 2*x0 % in linear zone
    n = 1e3; % barrier constant
    if x>=x0
        a = exp(n*(x-x0));
    else
        a = 2-exp(-n*x);
    end
    dmrur = a*dmrur;
end

B = [0, 0;
    0, 0;
    0, 0;
    1.0/J*torque_drive_c/omega_comp, 0;
    0, dmrur];

% Linearized C matrix
C = [0, 1, 0, 0, 0;
    p2/(SD_c(1)*p1^2), -1/(SD_c(1)*p1), 1, 0, 0];

C(2,:) = 100*C(2,:); % SD output is multiplied by factor 100



end


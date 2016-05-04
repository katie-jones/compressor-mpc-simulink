function [f,m_out] = get_comp_deriv(x, u, flag)
%#eml
% Inputs
torque_drive = u(1);  
Inflow_opening = u(2);
Outflow_opening = u(3); 
Recycle_opening = u(4);
In_pres = u(5); 
Out_pres = u(6);
q_in = -1;
if length(u)>6
 q_in = u(7);
end



% States
p1 = x(1);
p2 = x(2);
m_comp = x(3);
omega_comp = x(4);
m_rec = x(5);

[J,tauRecycle,A,C,m_in_c,m_rec_ss_c,D2,m_out_c,T_ss_c,~,torque_drive_c] = comp_coeffs();
[SpeedSound,~,~,VolumeT1,VolumeT2,AdivL] = const_flow();



torque_drive = torque_drive * torque_drive_c / (2 * pi * 50);
% if omega_comp > 2*pi*50
        torque_drive = torque_drive * (2 * pi * 50 / omega_comp);
% end

% Parameters
% NEW: Parameters for Air
% cv = 20.0;
% cp = 29.0;
% gamma = cp/cv; % 1.4028
% R_air = 8.31 / 28.97 *1000;
% T0 = 25;





% Algebraic equations

if In_pres > 0
dp_sqrt = sqrt(abs(In_pres*100 - p1*100)) * sign(In_pres*100 - p1*100);
M3 = [dp_sqrt*Inflow_opening^3 dp_sqrt*Inflow_opening^2 dp_sqrt*Inflow_opening dp_sqrt ...
    Inflow_opening^3 Inflow_opening^2 Inflow_opening 1]';

    m_in = C * M3 + m_in_c; % Inflow valve
else
    m_in = q_in;
end


% If using for simulation, have dead zone on recycle opening
% If using for MPC, remove dead zone (better results)
v = m_rec_ss_c(1) * sqrt(p2*1e5 - p1*1e5) * Recycle_opening  + m_rec_ss_c(2);
x0 = 1e-2;
if flag==1 || Recycle_opening >= 2*x0
    m_rec_ss = v * ~(Recycle_opening<x0);
else
    n = 1e2; % barrier constant
    if Recycle_opening>=x0
        a = 0.1+0.9*exp(n*(Recycle_opening-x0));
    else
        a = 2-0.9*exp(-n*Recycle_opening);
    end
    m_rec_ss = a * v;
end



dp_sqrt2 = sqrt(abs(p2*100 - Out_pres*100)) * sign(p2*100 - Out_pres*100);      

M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
       Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';

m_out = D2 * M5 + m_out_c;


% compressor pressure ratio

M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
    omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
    m_comp^3 m_comp^2 m_comp 1]';
p_ratio = A * M;


T_ss_model = T_ss_c(1) + T_ss_c(2) * m_comp;
T_ss_model = T_ss_model + T_ss_c(3);


f = zeros(5,1);
f(1) = SpeedSound * SpeedSound / VolumeT1 * (m_in + m_rec - m_comp) * 1e-5; % p1
f(2) = SpeedSound * SpeedSound / VolumeT2 * (m_comp - m_rec - m_out) * 1e-5;  % p2
f(3) = AdivL * (p_ratio * p1*1e5 - p2*1e5 ); % m_c
% sys(4) = tauComp * (dp - p_ratio); % phi
f(4) = 1.0 / J * (torque_drive - T_ss_model); % omega
f(5) = tauRecycle * (m_rec_ss - m_rec); % rec flow
% sys(7) = tauIn * (m_in_ss - m_in); % inflow
% sys(8) = tauOut * (m_out_ss - m_out); % outflow


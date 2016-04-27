function [sys,x0,str,tss]=comp_rig_model(t,x,u,flag,Param,X_ss)

switch flag,

case 0,	% Initialize the states and sizes
   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
    % ****************
  	%  Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states

   sys = mdlDerivatives(t,x,u,Param);

otherwise,
   sys = [];
end

% ******************************************
% Sub-routines or Functions
% ******************************************

% ******************************************
% Initialization
% ******************************************

function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss)



% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 5;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 10;     % outputs of model 
sizes.NumInputs      = 7;     % inputs of model
sizes.DirFeedthrough = 1;     % System is causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the states 


str = [];	                  % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

% ******************************************
%  Outputs
% ******************************************

function [sys] = mdlOutputs(t,x,u,Param)


% Inputs
% torque_drive = u(1);  
Inflow_opening = u(2);
Outflow_opening = u(3); 
% Recycle_opening = u(4);
% In_pres = u(5); 
Out_pres = u(6);
use_qin = u(7);

% States
p1 = x(1);%
p2 = x(2);%
m_comp = x(3);%
omega_comp = x(4);%
m_rec = x(5);%

[~,~,A,C,m_in_c,~,D2,m_out_c,T_ss_c,SD_c] = comp_coeffs();

% NEW: Parameters for Air
% cv = 20.0;
% cp = 29.0;
% gamma = cp/cv; % ca. 1.4028
% R_air = 8.31 / 28.97 * 1000;
% T_in = u(5);
% T0 = 25;


M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
    omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
    m_comp^3 m_comp^2 m_comp 1]';
dp = A * M;

% Map variation due to temperature:
% dp = ((1-T0/T_in) + T0/T_in * dp^((gamma-1)/gamma))^(gamma/(gamma-1));

% if dp < 1
%     dp = 1;
% end

if flag==0 % if mass flow in not given
    In_pres = u(5);
    dp_sqrt = sqrt(abs(In_pres*100 - p1*100)) * sign(In_pres*100 - p1*100);
    M3 = [dp_sqrt*Inflow_opening^3 dp_sqrt*Inflow_opening^2 dp_sqrt*Inflow_opening dp_sqrt ...
        Inflow_opening^3 Inflow_opening^2 Inflow_opening 1]';


    m_in = C * M3 + m_in_c; % Inflow valve
else
    m_in = u(5);
end


dp_sqrt2 = sqrt(abs(p2*100 - Out_pres*100)) * sign(p2*100 - Out_pres*100);      

M5 = [dp_sqrt2*Outflow_opening^3, dp_sqrt2*Outflow_opening^2, dp_sqrt2*Outflow_opening, dp_sqrt2, ...
       Outflow_opening^3, Outflow_opening^2, Outflow_opening, 1]';

m_out = D2 * M5 + m_out_c;

% output valve 
% m_out = Valve_out_gain*Outflow_opening*sqrt(abs(p2 - Out_pres*1e3))*sign(p2 - Out_pres*1e3); % Outflow valve
% m_out = max(m_out+0.05, 0);


T_ss_model = (T_ss_c(1) +  T_ss_c(2) * m_comp);% .* ~(m_comp < 1e-2);
% T_ss_model = signal_saturation(T_ss_model,0,40);
T_ss_model = T_ss_model + T_ss_c(3);


SD = -((p2*1e5) / (p1*1e5) / SD_c(1) - SD_c(2) / SD_c(1) - m_comp);


if isnan(dp) || isinf(dp) 
    keyboard;
end

sys(1) = p1; % p1
sys(2) = p2; % p2
sys(3) = m_comp; % m_comp 
sys(4) = dp; % p_ratio
sys(5) = omega_comp; % omega_comp
sys(6) = m_in; % 
sys(7) = m_out; % 
sys(8) = m_rec; % 
sys(9) = T_ss_model;
sys(10) = SD;

% keyboard;




% ******************************************
% Derivatives
% ******************************************


function sys = mdlDerivatives(t,x,u,Param)
sys = get_comp_deriv(x,u(1:6),u(7));


% Inputs
% torque_drive = u(1);  
% Inflow_opening = u(2);
% Outflow_opening = u(3); 
% Recycle_opening = u(4);
% dummy = u(5); 
% 
% 
% % States
% p1 = x(1);
% p2 = x(2);
% m_comp = x(3);
% omega_comp = x(4);
% m_rec = x(5);
% 
% [J,tauRecycle,A,C,m_in_c,m_rec_ss_c,D2,m_out_c,T_ss_c,~,torque_drive_c] = comp_coeffs();
% [SpeedSound,In_pres,Out_pres,VolumeT1,VolumeT2,AdivL] = const_flow();
% 
% torque_drive = torque_drive * torque_drive_c / (2 * pi * 50);
% % if omega_comp > 2*pi*50
%         torque_drive = torque_drive * (2 * pi * 50 / omega_comp);
% % end
% 
% % Parameters
% % NEW: Parameters for Air
% % cv = 20.0;
% % cp = 29.0;
% % gamma = cp/cv; % 1.4028
% % R_air = 8.31 / 28.97 *1000;
% % T0 = 25;
% 
% 
% dp_sqrt = sqrt(abs(In_pres*100 - p1*100)) * sign(In_pres*100 - p1*100);
% M3 = [dp_sqrt*Inflow_opening^3 dp_sqrt*Inflow_opening^2 dp_sqrt*Inflow_opening dp_sqrt ...
%     Inflow_opening^3 Inflow_opening^2 Inflow_opening 1]';
% 
% 
% m_in = C * M3 + m_in_c; % Inflow valve
% 
% m_rec_ss = ((m_rec_ss_c(1)) * sqrt(p2*1e5 - p1*1e5) * Recycle_opening  + m_rec_ss_c(2) ) * ~(Recycle_opening<1e-2);
% 
% 
% dp_sqrt2 = sqrt(abs(p2*100 - Out_pres*100)) * sign(p2*100 - Out_pres*100);      
% 
% M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
%        Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';
% 
% m_out = D2 * M5 + m_out_c;
% 
% % compressor pressure ratio
% 
% M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
%     omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
%     m_comp^3 m_comp^2 m_comp 1]';
% p_ratio = A * M;
% 
% % Map variation due to temperature:
% % dp = ((1-T0/T_in) + T0/T_in * dp^((gamma-1)/gamma))^(gamma/(gamma-1));
% 
% % if dp < 1
% %     dp = 1;
% % end
% 
% T_ss_model = (T_ss_c(1) +  T_ss_c(2) * m_comp); % .* ~(m_comp < 1e-2);
% % T_ss_model = signal_saturation(0,40);
% T_ss_model = T_ss_model + T_ss_c(3);
% 
% % if omega_comp > 2*pi*50
% %         T_ss_model = T_ss_model * (2 * pi * 50 / omega_comp);
% % end
% 
% sys(1) = SpeedSound * SpeedSound / VolumeT1 * (m_in + m_rec - m_comp) * 1e-5; % p1
% sys(2) = SpeedSound * SpeedSound / VolumeT2 * (m_comp - m_rec - m_out) * 1e-5;  % p2
% sys(3) = AdivL * (p_ratio * p1*1e5 - p2*1e5 ); % m_c
% % sys(4) = tauComp * (dp - p_ratio); % phi
% sys(4) = 1.0 / J * (torque_drive - T_ss_model); % omega
% sys(5) = tauRecycle * (m_rec_ss - m_rec); % rec flow
% % sys(7) = tauIn * (m_in_ss - m_in); % inflow
% % sys(8) = tauOut * (m_out_ss - m_out); % outflow
% 
% % NEW: Differential equation for T_mix
% % sys(9) = 1 / ( cv * p1 * VolumeT1 / (R_air * T_mix)) * (  cp * m_in * T_in + cp * m_rec * T_out -  cp * m_comp * T_mix ...
% %     - cv * T_mix * (m_in + m_rec - m_comp));
% % if abs(sys(9)) > 100
% %     sys(9) = 100;
% % end
% % keyboard
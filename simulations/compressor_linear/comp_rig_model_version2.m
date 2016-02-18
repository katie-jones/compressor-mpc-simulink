function [sys,x0,str,tss]=comp_rig_model_version2(t,x,u,flag,Param,X_ss)

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
sizes.NumInputs      = 5;     % inputs of model
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
torque_drive = u(1);  
Inflow_opening = u(2);
Outflow_opening = u(3); 
Recycle_opening = u(4);
dummy = u(5); 

% States
p1 = x(1);%
p2 = x(2);%
m_comp = x(3);%
omega_comp = x(4);%
m_rec = x(5);%

torque_drive = torque_drive * 15000 / (2 * pi * 50);
% if omega_comp > 2*pi*50
        torque_drive = torque_drive * (2 * pi * 50 / omega_comp);
% end

% NEW: Parameters for Air
% cv = 20.0;
% cp = 29.0;
% gamma = cp/cv; % ca. 1.4028
% R_air = 8.31 / 28.97 * 1000;
% T_in = u(5);
% T0 = 25;

SpeedSound = 340;     % speed of sound
In_pres    = 1.0;  % input pressure
Out_pres   = 1.0; % output pressure



%%% coefficients of the compressor map 
A = [0.000299749505193654     -0.000171254191089237      3.57321648097597e-05    -9.1783572200945e-07 ...
    -0.252701086129365         0.136885752773673         -0.02642368327081        0.00161012740365743 ...
    54.8046725371143          -29.9550791497765          5.27827499839098         0.693826282579158];

%%% coefficients of the torque map 
B = [0.0230082553748617       -0.00142888713858701     36.1012418352137          3.34780616372668];

M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
    omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
    m_comp^3 m_comp^2 m_comp 1]';
dp = A * M;

% Map variation due to temperature:
% dp = ((1-T0/T_in) + T0/T_in * dp^((gamma-1)/gamma))^(gamma/(gamma-1));

% if dp < 1
%     dp = 1;
% end


% compressor torque
% M2 = [m_comp*omega_comp omega_comp m_comp 1]';
% torque_comp =  B * M2;
% torque_comp = signal_saturation(torque_comp,0,40);

C = [-0.423884232813775         0.626400271518973       -0.0995040168384753      0.0201535563630318 ...
     -0.490814924104294         0.843580880467905        -0.423103455111209      0.0386841406482887];

dp_sqrt = sqrt(abs(In_pres*100 - p1*100)) * sign(In_pres*100 - p1*100);
M3 = [dp_sqrt*Inflow_opening^3 dp_sqrt*Inflow_opening^2 dp_sqrt*Inflow_opening dp_sqrt ...
    Inflow_opening^3 Inflow_opening^2 Inflow_opening 1]';


m_in = C * M3 + 0.0051; % Inflow valve
% m_rec_ss = Valve_rec_gain*Recycle_opening*sqrt(abs(p2 - p1))*sign(p2 - p1); % Recycle valve
% 
% D = [         0.021182190468223
%                          0
%        -0.0372607070498975
%          0.015018502625883
%         -0.163754411108535
%                          0
%          0.395432975110223
%         -0.129295806042816
%          0.036946153978411
%                          0
%        -0.0787044830112647
%         0.0458965683496126]';
% dp_sqrt2 = sqrt(abs(p2/1000 - Out_pres)) .* sign(p2/1000 - Out_pres);     
% M4 = [dp_sqrt2.^2.*Outflow_opening.^3 dp_sqrt2.^2.*Outflow_opening.^2 dp_sqrt2.^2.*Outflow_opening dp_sqrt2.^2 ...
%       dp_sqrt2.*Outflow_opening.^3 dp_sqrt2.*Outflow_opening.^2 dp_sqrt2.*Outflow_opening dp_sqrt2 ...
%       Outflow_opening.^3 Outflow_opening.^2 Outflow_opening 1]';
% 
% m_out = D * M4;
D2 = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];
dp_sqrt2 = sqrt(abs(p2*100 - Out_pres*100)) * sign(p2*100 - Out_pres*100);      
% M4 = [dp_sqrt2.^2.*u_out.^3 dp_sqrt2.^2.*u_out.^2 dp_sqrt2.^2.*u_out dp_sqrt2.^2 ...
%       dp_sqrt2.*u_out.^3 dp_sqrt2.*u_out.^2 dp_sqrt2.*u_out dp_sqrt2 ...
%       u_out.^3 u_out.^2 u_out ones(length(u_out),1)]';
M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
       Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';
% m_out = D * M4 + 0.05;
m_out = D2 * M5 + 0.0170;

% output valve 
% m_out = Valve_out_gain*Outflow_opening*sqrt(abs(p2 - Out_pres*1e3))*sign(p2 - Out_pres*1e3); % Outflow valve
% m_out = max(m_out+0.05, 0);


T_ss_model = (2.5543945754982 +  47.4222669576423 * m_comp);% .* ~(m_comp < 1e-2);
% T_ss_model = signal_saturation(T_ss_model,0,40);
T_ss_model = T_ss_model + 0.6218;


SD = -((p2*1e5) / (p1*1e5) / 5.55 - 0.66 / 5.55 - m_comp);


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


% Inputs
torque_drive = u(1);  
Inflow_opening = u(2);
Outflow_opening = u(3); 
Recycle_opening = u(4);
dummy = u(5); 


% States
p1 = x(1);
p2 = x(2);
m_comp = x(3);
omega_comp = x(4);
m_rec = x(5);

torque_drive = torque_drive * 15000 / (2 * pi * 50);
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

SpeedSound = 340;     % speed of sound
In_pres    = 1.0;  % input pressure
Out_pres   = 1.0; % output pressure

% Valve_in_gain = 1 / 270; % 420.0 / 3600 / sqrt(9789 / 11.82);   % opening gain
% Valve_rec_gain = 420.0 / 3600 / sqrt(9789 / 11.82); %0.01; % opening gain
% Valve_out_gain = 1 / 270; % 420.0 / 3600 / sqrt(9789 / 11.82); %1.0/sqrt(0.2e5); % opening gain

VolumeT1 = 2*pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 8.191; % volume of the inlet pipe, 6m
VolumeT2 =  pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 5.940; % volume, tank
AdivL =     pi * (0.08 / 2)^2 / 3 * 0.1;  % Duct cross section divided by length

%%% Compressor
J          = (0.4 + 0.2070) *0.4; % 0.8;  % shaft inertia: compressor and rotor
% tauComp    = 1.0/0.001;     % time constant of the pressure ratio
tauRecycle =   1/0.5 * 1;     % time constant of the recycle valve




%%% coefficients of the compressor map 
A = [0.000299749505193654     -0.000171254191089237      3.57321648097597e-05    -9.1783572200945e-07 ...
    -0.252701086129365         0.136885752773673         -0.02642368327081        0.00161012740365743 ...
    54.8046725371143          -29.9550791497765          5.27827499839098         0.693826282579158];

%%% coefficients of the torque map 
B = [0.0230082553748617       -0.00142888713858701     36.1012418352137          3.34780616372668];


% Algebraic equations

% valves
C = [-0.423884232813775         0.626400271518973       -0.0995040168384753      0.0201535563630318 ...
     -0.490814924104294         0.843580880467905        -0.423103455111209      0.0386841406482887];

dp_sqrt = sqrt(abs(In_pres*100 - p1*100)) * sign(In_pres*100 - p1*100);
M3 = [dp_sqrt*Inflow_opening^3 dp_sqrt*Inflow_opening^2 dp_sqrt*Inflow_opening dp_sqrt ...
    Inflow_opening^3 Inflow_opening^2 Inflow_opening 1]';


m_in = C * M3 + 0.0051; % Inflow valve
% m_rec_ss = Valve_rec_gain*Recycle_opening*sqrt(abs(p2 - p1))*sign(p2 - p1); % Recycle valve
% if Recycle_opening > 1e-2
%     REC = [-0.70259 0.69713 -0.02843 0.013383 3.6089 -4.2137 1.3645 -0.11049];
%     dp_sqrt8 = sqrt(abs(p2/1000 - p1/1000));
%     M3 = [dp_sqrt8.*Recycle_opening.^3 dp_sqrt8.*Recycle_opening.^2 dp_sqrt8.*Recycle_opening dp_sqrt8 ...
%         Recycle_opening.^3 Recycle_opening.^2 Recycle_opening ones(length(Recycle_opening),1)]';
%     m_rec_ss = REC * M3 + 0.05;
% else
%     m_rec_ss = 0;
% end
m_rec_ss = ((0.0047) * sqrt(p2*1e5 - p1*1e5) * Recycle_opening  + 0.0263 ) * ~(Recycle_opening<1e-2);


% D = [    0.021182190468223
%                          0
%        -0.0372607070498975
%          0.015018502625883
%         -0.163754411108535
%                          0
%          0.395432975110223
%         -0.129295806042816
%          0.036946153978411
%                          0
%        -0.0787044830112647
%         0.04589656834961267]';
% dp_sqrt2 = sqrt(abs(p2/1000 - Out_pres)) .* sign(p2/1000 - Out_pres);     
% M4 = [dp_sqrt2.^2.*Outflow_opening.^3 dp_sqrt2.^2.*Outflow_opening.^2 dp_sqrt2.^2.*Outflow_opening dp_sqrt2.^2 ...
%       dp_sqrt2.*Outflow_opening.^3 dp_sqrt2.*Outflow_opening.^2 dp_sqrt2.*Outflow_opening dp_sqrt2 ...
%       Outflow_opening.^3 Outflow_opening.^2 Outflow_opening 1]';
% 
% m_out = D * M4;
% % output valve 
% % m_out = Valve_out_gain*Outflow_opening*sqrt(abs(p2 - Out_pres*1e3))*sign(p2 - Out_pres*1e3); % Outflow valve
% m_out = max(m_out+0.05, 0);
D2 = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];
dp_sqrt2 = sqrt(abs(p2*100 - Out_pres*100)) * sign(p2*100 - Out_pres*100);      
% M4 = [dp_sqrt2.^2.*u_out.^3 dp_sqrt2.^2.*u_out.^2 dp_sqrt2.^2.*u_out dp_sqrt2.^2 ...
%       dp_sqrt2.*u_out.^3 dp_sqrt2.*u_out.^2 dp_sqrt2.*u_out dp_sqrt2 ...
%       u_out.^3 u_out.^2 u_out ones(length(u_out),1)]';
M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
       Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';
% m_out = D * M4 + 0.05;
m_out = D2 * M5 + 0.0170;

% compressor pressure ratio

M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
    omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
    m_comp^3 m_comp^2 m_comp 1]';
p_ratio = A * M;

% Map variation due to temperature:
% dp = ((1-T0/T_in) + T0/T_in * dp^((gamma-1)/gamma))^(gamma/(gamma-1));

% if dp < 1
%     dp = 1;
% end

% compressor torque
M2 = [m_comp*omega_comp omega_comp m_comp 1]';
% torque_comp =  B * M2;
T_ss_model = (2.5543945754982 +  47.4222669576423 * m_comp); % .* ~(m_comp < 1e-2);
% T_ss_model = signal_saturation(0,40);
T_ss_model = T_ss_model + 0.6218;

% if omega_comp > 2*pi*50
%         T_ss_model = T_ss_model * (2 * pi * 50 / omega_comp);
% end

sys(1) = SpeedSound * SpeedSound / VolumeT1 * (m_in + m_rec - m_comp) * 1e-5; % p1
sys(2) = SpeedSound * SpeedSound / VolumeT2 * (m_comp - m_rec - m_out) * 1e-5;  % p2
sys(3) = AdivL * (p_ratio * p1*1e5 - p2*1e5 ); % m_c
% sys(4) = tauComp * (dp - p_ratio); % phi
sys(4) = 1.0 / J * (torque_drive - T_ss_model); % omega
sys(5) = tauRecycle * (m_rec_ss - m_rec); % rec flow
% sys(7) = tauIn * (m_in_ss - m_in); % inflow
% sys(8) = tauOut * (m_out_ss - m_out); % outflow

% NEW: Differential equation for T_mix
% sys(9) = 1 / ( cv * p1 * VolumeT1 / (R_air * T_mix)) * (  cp * m_in * T_in + cp * m_rec * T_out -  cp * m_comp * T_mix ...
%     - cv * T_mix * (m_in + m_rec - m_comp));
% if abs(sys(9)) > 100
%     sys(9) = 100;
% end
% keyboard
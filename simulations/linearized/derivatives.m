%% dp1dot/dp1
clear all
u = zeros(5);
x = zeros(5);
dp = 0.02;

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

% Calculate p1dot varying p1
for i=1:100
    p1 = -1+i*dp;

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

% p1dot
pdot1(i) = sys(1);

%dp1dot/dp1
dx1x1(i) = SpeedSound * SpeedSound / VolumeT1 * 1e-5 * (-1/2*sign(In_pres-p1)*100/sqrt(abs(In_pres*100-p1*100))) * (C(1)*Inflow_opening^3 + C(2)*Inflow_opening^2 + C(3)*Inflow_opening + C(4));
end

% estimate of dp1dot/dp1 using forward euler
dx1x1_est = (pdot1(2:end)-pdot1(1:end-1))/dp;
figure; hold on
plot(dx1x1)
plot(dx1x1_est,'-g')


%% dp1dot/dqc and dp1dot/dqr
dx1x3 = -SpeedSound * SpeedSound / VolumeT1 * 1e-5;
dx1x5 = SpeedSound * SpeedSound / VolumeT1 * 1e-5;

%% dp2dot/dp2

% Calculate p2dot varying p2
for i=1:100
    p2 = -1+i*dp;

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

% p2dot
pdot2(i) = sys(2);

%dp2dot/dp2
dx2x2(i) = SpeedSound * SpeedSound / VolumeT2 * 1e-5 * (1/2*sign(p2-Out_pres)*100/sqrt(abs(p2*100-Out_pres*100))) * (D2(1)*Outflow_opening^3 + D2(2)*Outflow_opening^2 + D2(3)*Outflow_opening + D2(4));
end

% estimate of dp1dot/dp1 using forward euler
dx2x2_est = (pdot2(2:end)-pdot2(1:end-1))/dp;
figure; hold on
plot(dx2x2)
plot(dx2x2_est,'-g')


%% dp2dot/dqc and dqr

dx2x5 = -SpeedSound * SpeedSound / VolumeT2 * 1e-5;
dx2x3 = SpeedSound * SpeedSound / VolumeT2 * 1e-5;

%% dqcdot/dqc
for i=1:100
    m_comp = -1+i*dp;

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

% p2dot
qcdot(i) = sys(3);

%dp3dot/dp3
dM_dmcomp = [3*omega_comp^2*m_comp^2, 2*omega_comp^2*m_comp, omega_comp^2, 0,...
            3*omega_comp*m_comp^2, 2*omega_comp*m_comp, omega_comp, 0,...
            3*m_comp^2, 2*m_comp, 1, 0]';
dx3x3(i) = AdivL * (p1*1e5)*A*dM_dmcomp;

end

dx3x3_est = (qcdot(2:end)-qcdot(1:end-1))/dp;

figure; hold on;
plot(dx3x3);
plot(dx3x3_est,'-g')


%% dqcdot/dp1 and /dp2
dx3x1 = AdivL*(p_ratio*1e5);
dx3x2 = -AdivL*1e5;

%% dqcdot/domegac
for i=1:100
    omega_comp = -1+i*dp;
    m_comp = 1;

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

% p2dot
qcdot(i) = sys(3);

%dp3dot/dp3
dM_dwcomp = [2*omega_comp*m_comp^3, 2*omega_comp*m_comp^2, 2*omega_comp*m_comp, 2*omega_comp,...
            m_comp^3, m_comp^2, m_comp, 1,...
            0, 0, 0, 0]';
dx3x4(i) = AdivL * (p1*1e5)*A*dM_dwcomp;

end

dx3x4_est = (qcdot(2:end)-qcdot(1:end-1))/dp;

figure; hold on;
plot(dx3x4);
plot(dx3x4_est,'-g')


%% domegacdot/dTd and /dqc
dx4du1 = 1.0/J;
dx4dx3 = -1.0/J * 47.4222669576423;

%% dqrdot/dqr and /dropen
dx5x5 = -tauRecycle;
dx5u2 = tauRecycle*0.0047 * sqrt(p2*1e5 - p1*1e5);

%% dqrdot/dp1 and /dp2

for i=1:100
    p2 = 1+i*dp;
    Recycle_opening = 0.1;
    p1 = 0;

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

% p2dot
qrdot(i) = sys(5);

dx5x1(i) = -tauRecycle * (0.0047 * 1/2 * Recycle_opening / sqrt(p2*1e5 - p1*1e5) * 1e5);
dx5x2(i) = tauRecycle * (0.0047 * 1/2 * Recycle_opening / sqrt(p2*1e5 - p1*1e5) * 1e5);
end

dx5x1_est = (qrdot(2:end)-qrdot(1:end-1))/dp;

figure; hold on;
plot(dx5x2);
plot(dx5x1_est,'-g')

%% dSD/dp2 and /dp1 and /dqc
dy2x1 = p2/(5.55*p1^2);
dy2x2 = -1/(5.55*p1);
dy2x3 = 1;

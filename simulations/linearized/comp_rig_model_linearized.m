function [sys,x0,str,tss]=comp_rig_model_linearized(t,x,u,flag,Param,X_ss)
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
% global x_old u_old
% 
% if ~(exist('x_old','var'))
%     x_old = x;
%     u_old = zeros(5,1);
% end

u_old = [    
    0.3040
%     0
    0.4050
    0.3930
         0
         0];
x_old = [ 0.8980
    1.1260
    0.1500
  439.5000
         0];

f = get_comp_deriv(x_old,u_old(1:4))';
[A,B,C] = get_linearized_matrices(x_old,u_old([1,4],1));

du = (u-u_old);
xdot = A*(x-x_old) + B*(du([1,4],1)) + f;
sys = xdot;

x_old = x;
u_old = u;



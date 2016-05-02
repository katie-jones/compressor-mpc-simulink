function [sys,x0,str,tss]=comp_rig_parallel_linearized(t,x,u,flag,Param,X_ss)
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
sizes.NumContStates  = 11;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 21;     % outputs of model 
sizes.NumInputs      = 12;     % inputs of model
sizes.DirFeedthrough = 1;     % System is causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the states 


str = [];	                  % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].

end

% ******************************************
%  Outputs
% ******************************************

function [sys] = mdlOutputs(t,x,u,Param)
[~,xsize,~,usize] = const_sim;

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);
PD = x(end);

u1(end) = PD;
u2(end) = PD;

sys_comp1 = compOutputs(x1,u1);
sys_comp2 = compOutputs(x2,u2);

sys = [sys_comp1(:); sys_comp2(:); PD];

end

function sys = compOutputs(x,u)

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

[J,tauRecycle,A,C,m_in_c,m_rec_ss_c,D2,m_out_c,T_ss_c,SD_c,torque_drive_c] = comp_coeffs();
[SpeedSound,In_pres,Out_pres,VolumeT1,VolumeT2,AdivL] = const_flow();

% for parallel simulation with non-consant output pressure
if dummy > 0
    Out_pres = dummy;
end

torque_drive = torque_drive * torque_drive_c / (2 * pi * 50);
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


M = [omega_comp^2*m_comp^3 omega_comp^2*m_comp^2 omega_comp^2*m_comp omega_comp^2 ...
    omega_comp*m_comp^3 omega_comp*m_comp^2 omega_comp*m_comp omega_comp ...
    m_comp^3 m_comp^2 m_comp 1]';
dp = A * M;

% Map variation due to temperature:
% dp = ((1-T0/T_in) + T0/T_in * dp^((gamma-1)/gamma))^(gamma/(gamma-1));

% if dp < 1
%     dp = 1;
% end


dp_sqrt = sqrt(abs(In_pres*100 - p1*100)) * sign(In_pres*100 - p1*100);
M3 = [dp_sqrt*Inflow_opening^3 dp_sqrt*Inflow_opening^2 dp_sqrt*Inflow_opening dp_sqrt ...
    Inflow_opening^3 Inflow_opening^2 Inflow_opening 1]';


m_in = C * M3 + m_in_c; % Inflow valve


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

end


% ******************************************
% Derivatives
% ******************************************


function sys = mdlDerivatives(t,x,u,Param)
global x_old u_old

[~,xsize,~,usize,~,uoff1,uoff2,ud] = const_sim;


if (isempty(x_old))
    x_init_lin = [0.899; 1.126; 0.15; 440; 0];
    x_old = [x_init_lin; x_init_lin; 1.08];
    u_old = [uoff1; uoff2; 0.3; ud];
end
    

u1 = u_old(1:usize);
u2 = u_old(usize+1:2*usize);

x1 = x_old(1:xsize);
x2 = x_old(xsize+1:2*xsize);
PD = x_old(end);

u1(end) = PD;
u2(end) = PD;

ud = u(2*usize+1);
m_tank_in = u(end);

dx1 = get_comp_deriv(x1,u1,1);
dx2 = get_comp_deriv(x2,u2,1);
dtank = get_tank_deriv(PD,[m_tank_in; ud]);

du = u - u_old;

[A,B] = linearize_tank(x_old, [u1; u2; ud]);

sys = [dx1; dx2; dtank] + A*(x-x_old) + B*[du(1); du(4); du(1+usize); du(4+usize)];

u_old = u;
x_old = x;

end


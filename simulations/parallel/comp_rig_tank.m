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
sizes.NumContStates  = 1;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 2;     % outputs of model 
sizes.NumInputs      = 2;     % inputs of model
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
global P_D

% States
P_D = x(1);

% Inputs
m_in = u(1);
Outflow_opening = u(2);

SpeedSound = 340;     % speed of sound
Out_pres = 1.0; % atmospheric pressure


D2 = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];
dp_sqrt2 = sqrt(abs(P_D*100 - Out_pres*100)) * sign(P_D*100 - Out_pres*100);      
% M4 = [dp_sqrt2.^2.*u_out.^3 dp_sqrt2.^2.*u_out.^2 dp_sqrt2.^2.*u_out dp_sqrt2.^2 ...
%       dp_sqrt2.*u_out.^3 dp_sqrt2.*u_out.^2 dp_sqrt2.*u_out dp_sqrt2 ...
%       u_out.^3 u_out.^2 u_out ones(length(u_out),1)]';
M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
       Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';
% m_out = D * M4 + 0.05;
m_out = D2 * M5 + 0.0170;

sys(1) = m_out; % p1
sys(2) = P_D; % p2


% ******************************************
% Derivatives
% ******************************************


function sys = mdlDerivatives(t,x,u,Param)
global P_D

% States
P_D = x(1);

% Inputs
m_in = u(1);
Outflow_opening = u(2);

SpeedSound = 340;     % speed of sound
Out_pres = 1.0; % atmospheric pressure

VolumeT =  pi * (0.60 / 2)^2 * 2 + pi * (0.08 / 2)^2 * 5.940; % volume, tank

D2 = [-0.0083454 -0.0094965 0.16826 -0.032215 -0.61199 0.94175 -0.48522 0.10369];
dp_sqrt2 = sqrt(abs(P_D*100 - Out_pres*100)) * sign(P_D*100 - Out_pres*100);      
% M4 = [dp_sqrt2.^2.*u_out.^3 dp_sqrt2.^2.*u_out.^2 dp_sqrt2.^2.*u_out dp_sqrt2.^2 ...
%       dp_sqrt2.*u_out.^3 dp_sqrt2.*u_out.^2 dp_sqrt2.*u_out dp_sqrt2 ...
%       u_out.^3 u_out.^2 u_out ones(length(u_out),1)]';
M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
       Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';
% m_out = D * M4 + 0.05;
m_out = D2 * M5 + 0.0170;

sys(1) = SpeedSound * SpeedSound / VolumeT * (m_in - m_out) * 1e-5; % p1


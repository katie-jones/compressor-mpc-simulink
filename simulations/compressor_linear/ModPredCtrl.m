function [sys,x0,str,tss]=ModPredCtrl(t,x,u,flag,Param,X_ss)

global xinit H A B C M ysize xsize usize p LB UB Ga Gb Gc Ain b upast m

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
function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);

global xinit H A B C M ysize xsize usize p LB UB Ga Gb Gc Ain b upast m

% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 0;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 2;     % outputs of model
sizes.NumInputs      = 4;     % inputs of model
sizes.DirFeedthrough = 1;     % system is causal, there is explicit dependence of outputs on inputs
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % initialize the discrete states
str = [];	                  % set str to an empty matrix
tss = [0.05,0];	              % sample time: [period, offset]

% ******************************************
%  Outputs
% ******************************************
function sys = mdlOutputs(t,x,u,Param);

global xinit H A B C M ysize xsize usize p LB UB Ga Gb Gc Ain b upast m


% get the measurements for outputs and the output reference

% keyboard;

y=u(1:2);
yref=u(3:4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% estimate the state information


xinit  = xinit + M*(y - C*xinit);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upper and lower bounds wrt past input

for i=1:m
    Upast(1+(i-1)*usize:i*usize,1)=upast;
end

LBa=LB-Upast;
UBa=UB-Upast;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate the output reference vector


for i=1:p
    R(1+(i-1)*ysize:i*ysize,1)=yref;
end

% generate the gradient vector

f = xinit'*Gb+Upast'*Gc-R'*Ga;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% solve the QP
% keyboard;
% usol=quadprog(H,f,Ain,b,[],[],LBa,UBa);
% (H,g,A,lb,ub,lbA,ubA)
usol =  qpOASES(H,f',Ain,LBa,UBa,[],b);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% extract the control moves for next time step


dunext=usol(1:usize,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initial conditions are updated

% keyboard
xinit= A*xinit + B*(dunext+upast);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control moves are sent to the plant

upast=dunext+upast;

sys=upast;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
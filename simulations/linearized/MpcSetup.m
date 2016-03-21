% Tuning parameters for MPC controller
clear all
% addpath('..\qpoases\interfaces\matlab\')
addpath('../call_qpoases_m')
addpath('../call_qpoases_m/qpoases3/interfaces/matlab')
addpath('../common')

% global xinit A B C M ysize xsize usize p LB UB Ain b upast m

load sys.mat

[xsizea,ysize,dsize,usize,n_delay,xsize,Ts,p,m,UWT,YWT] = mpc_constants();



%% Define augmented system
% delay of n_delay time steps in 2nd component of u
% constant output disturbance

Adelay1 = [zeros(n_delay(1)-1,1), eye(n_delay(1)-1); zeros(1,n_delay(1))]; % delay component of A
Adelay2 = [zeros(n_delay(2)-1,1), eye(n_delay(2)-1); zeros(1,n_delay(2))]; % delay component of A
Adist = eye(dsize); % disturbance component of A

if (n_delay(1)==0)
    Aaug = [A, B(:,2), zeros(xsizea,n_delay(2)-1), zeros(xsizea,dsize);
        zeros(n_delay(2),xsizea), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsizea), zeros(dsize,n_delay(2)), Adist];

    Baug = [B(:,1), zeros(xsizea,1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1];
        zeros(dsize,2)];
else
    Aaug = [A, B(:,1), zeros(xsizea,n_delay(1)-1), B(:,2), zeros(xsizea,n_delay(2)-1), zeros(xsizea,dsize);
        zeros(n_delay(1),xsizea), Adelay1, zeros(n_delay(1),n_delay(2)), zeros(n_delay(1),dsize);
        zeros(n_delay(2),xsizea), zeros(n_delay(2),n_delay(1)), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsizea), zeros(dsize,n_delay(1)), zeros(dsize,n_delay(2)), Adist];

    Baug = [zeros(xsizea,usize);
        [zeros(n_delay(1)-1,1); 1], zeros(n_delay(1),1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1); 1];
        zeros(dsize,usize)];
end

% disturbance to output matrix
Cd = eye(ysize,dsize); 

Caug = [C, zeros(ysize,sum(n_delay)), Cd];

Daug = 0;

sys = ss(Aaug,Baug,Caug,Daug,Ts);

%% Define correct global variables
A = Aaug;
B = Baug;
C = Caug;
orig_xsize = xsizea;
% xsize = length(Aaug);

%% Design observer
% expectation of output disturbance
Qn = eye(orig_xsize);
Rn = eye(dsize);

% state disturbance to state matrix
G = [zeros(xsize-dsize,orig_xsize);
    zeros(dsize,1), eye(dsize)*10, zeros(dsize,orig_xsize-dsize-1)];

Hkalman = [zeros(dsize,orig_xsize-dsize), eye(dsize)];

sys_kalman = ss(A,[B G], C, [D Hkalman], Ts);

[KEST,L,P,M,Z] = kalman(sys_kalman,Qn,Rn); % get observer values


%% Define upper/lower bounds
lb = [-0.3; 0];
ub = [0.3; 1];
% lb = [-0.1; 0];
% ub = [0.1; 1];
% lb = [0;0];
% ub = [0;1];
LB = repmat(lb,m,1);
UB = repmat(ub,m,1);

lb_rate = [0.1; 0.1];
ub_rate = [0.1; 1];
% lb_rate = [0.05;0.001];
% ub_rate = [0.05;0.5];
LBrate = repmat(lb_rate,m,1);
UBrate = repmat(ub_rate,m,1);

% A matrix for qpOASES
Ain = full(spdiags(ones(usize*m,1)*[1,-1],[-usize,0],usize*m,usize*m));

% replicate A to only use upper bound ubA
Ain = [Ain; -Ain];
b = [LBrate; UBrate];


%% initialize system state

% linearization point
x_init_lin = [0.898
      1.126
      0.15
      439.5
      0];

xinit = [x_init_lin;
        zeros(xsize-orig_xsize,1)];
upast = zeros(usize,1);
deltax = zeros(xsize,1);

disp('MPC Problem Formulated!');

generate_file_linearized;

disp('Embedded files generated.')

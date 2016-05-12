% Tuning parameters for MPC controller
clear all

addpath('../call_qpoases_m')
addpath('../call_qpoases_m/qpoases3/interfaces/matlab/')
addpath('../common')
addpath('../parallel_common')

%% Parameters for saving results
% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
results_folder = '../results';
results_label = 'Centralized';
results_fname = 'cen';
results_overwrite = 0;

if ~exist(results_folder,'dir')
    mkdir(results_folder)
end

saveplots = 1;

%% Constants
[Ts, xsize_comp, xsize, ~, ysize, uoff1, uoff2, ud] = const_sim();
[n_delay,dsize,ucontrolsize,p,m] = const_mpc();

xtotalsize = xsize + 2*sum(n_delay) + 2*dsize;
weights;

%% Initial state
P_D = 1.12;

x_init_lin = [0.916; 1.145; 0.152; 440; 0];

xinit = [x_init_lin; x_init_lin; P_D];

u_out = uoff1(3);
u_d = ud;

u_init = zeros(2*ucontrolsize,1);

[A,B,C] = get_qp_matrices(xinit, u_init, UWT, YWT);

xinit = [x_init_lin; x_init_lin; P_D];

D = zeros(ysize);

% Reference output
yss = [0.2175 0.2175 0 1.12]';
yref = [0.22 0.22 0 1.12]';

%% Design observer
% expectation of output disturbance
Qn = eye(xsize);
Rn = eye(2*dsize);

% state disturbance to state matrix
G = [zeros(xtotalsize-2*dsize,xsize);
    zeros(dsize,1), eye(dsize)*10, zeros(dsize, xsize_comp-dsize-1), zeros(dsize,xsize-xsize_comp);
    zeros(dsize,xsize_comp), zeros(dsize,1), eye(dsize)*10, zeros(dsize, xsize-xsize_comp-dsize-1)];

Hkalman = [zeros(dsize, xsize_comp-dsize), eye(dsize), zeros(dsize,xsize-xsize_comp);
    zeros(dsize,xsize_comp), zeros(dsize,xsize_comp-dsize), eye(dsize), zeros(dsize, xsize-2*xsize_comp)];

sys_kalman = ss(A, [B G], C, [D Hkalman], Ts);

[KEST, L, P, M, Z] = kalman(sys_kalman, Qn, Rn);

%% Define upper/lower bounds
lb = [-0.3; 0];
ub = [0.3; 1];
% lb = [-0.1; 0];
% ub = [0.1; 1];
% lb = [0;0];
% ub = [0;0];


lb_rate = [0.1; 0.1];
ub_rate = [0.1; 1];

LB = repmat(lb,2*m,1);
UB = repmat(ub,2*m,1);
LBrate = repmat(lb_rate,2*m,1);
UBrate = repmat(ub_rate,2*m,1);

usize = 2*ucontrolsize;
Ain = full(spdiags(ones(usize*m,1)*[1,-1],[-usize,0],usize*m,usize*m));

Ain = [Ain; -Ain];
b = [LBrate; UBrate];

%% Initialize system state
xinit = [xinit; zeros(xtotalsize-xsize,1)];
upast = u_init;
deltax = zeros(xtotalsize,1);

%% Define disturbances
for n_disturbance=1:4
[tdist,udist1,udist2,dist_dirname] = disturbances(n_disturbance);


disp('MPC Problem Formulated');

generate_file_linearized;

disp('Embedded files generated.');

sim('closedloop.mdl');
makeplots;
end
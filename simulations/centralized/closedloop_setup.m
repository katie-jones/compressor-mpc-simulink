% Tuning parameters for MPC controller
% clear all
addpath('../call_qpoases_m')
addpath('../call_qpoases_m/qpoases3/interfaces/matlab/')
addpath('../common')

%% Results variables
% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
saveplots = 0;

results_folder = '../results';
results_fname = 'cen';
results_overwrite = 0;

if ~exist(results_folder,'dir')
    mkdir(results_folder)
end

%% Constants
% Reference output
xinit = [  0.86733   1.0308  0.17648   394.96        0  0.99889   1.1872  0.17648   394.96        0 ]';

yss = [   1.0308   8.1258   1.1872   8.1258 ]';

yref = [1.03 8.12 1.1872 8.12];

[Ts, xsize_comp, xsize, ~, ysize, uoff1, uoff2] = const_sim();
[n_delay,dsize,ucontrolsize,p,m] = const_mpc();

xtotalsize = xsize + 2*sum(n_delay) + 2*dsize;

weights;


%% Initial state
% global P_D

u_out = uoff1(3);

u_init = zeros(2*ucontrolsize,1);

[A,B,C] = get_qp_matrices(xinit, u_init, UWT, YWT);


D = zeros(ysize,2*ucontrolsize);

%% Design observer
% expectation of output disturbance
Qn = eye(xsize);
Rn = eye(2*dsize);

% state disturbance to state matrix
G = [zeros(xtotalsize-2*dsize,xsize);
    zeros(dsize,1), eye(dsize)*10, zeros(dsize, xsize_comp-dsize-1), zeros(dsize,xsize-xsize_comp);
    zeros(dsize,xsize_comp), zeros(dsize,1), eye(dsize)*10, zeros(dsize, xsize-xsize_comp-dsize-1)];

Hkalman = [zeros(dsize, xsize_comp-dsize), eye(dsize), zeros(dsize,xsize-xsize_comp);
    zeros(dsize,xsize_comp), zeros(dsize,xsize_comp-dsize), eye(dsize), zeros(dsize, xsize-2*xsize_comp);
    zeros(ysize-2*dsize,xsize)];

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


disp('MPC Problem Formulated');

generate_file_linearized;

disp('Embedded files generated.');

%% Define disturbances
% Choose type of disturbance
% 1: output, 2: input, 3: asymmetric output, 4: asymmetric input, 5: big output
for n_disturbance=1
    [tdist,udist1,udist2,dist_dirname] = disturbances(n_disturbance);
%     sim('closedloop');
%     makeplots;
end

% Tuning parameters for MPC controller

addpath('../call_qpoases_m')
addpath('../call_qpoases_m/qpoases3/interfaces/matlab/')
addpath('../common')

[Ts, xsize_comp, xsize, ~, ysize, uoff1, uoff2] = const_sim();
[n_delay,dsize,ucontrolsize,p,m] = const_mpc();

xtotalsize = xsize + 2*sum(n_delay) + 2*dsize;


%% Initial state
u_init = zeros(2*ucontrolsize,1);
xinit = [xinit; zeros(xtotalsize-xsize,1)];
[A,B,C] = get_qp_matrices(xinit, u_init, zeros(ysize,1),UWT,YWT);

D = zeros(ysize);

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

LB = repmat(lb,m,1);
UB = repmat(ub,m,1);
LBrate = repmat(lb_rate,m,1);
UBrate = repmat(ub_rate,m,1);

usize = ucontrolsize;
Ain = full(spdiags(ones(usize*m,1)*[1,-1],[-usize,0],usize*m,usize*m));

Ain = [Ain; -Ain];
b = [LBrate; UBrate];

%% Initialize system state
upast = u_init;
deltax = zeros(xtotalsize,1);

% add disturbances
[tdist,udist1,udist2,dist_dirname] = disturbances(n_disturbance);

disp('MPC Problem Formulated');

generate_file_linearized;

disp('Embedded files generated.');


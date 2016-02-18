% Tuning parameters for MPC controller
% clear all
% addpath('..\qpoases\interfaces\matlab\')
addpath('call_qpoases_m')
addpath('call_qpoases_m/qpoases3/interfaces/matlab')

global xinit H A B C M ysize xsize usize p LB UB Ga Gb Gc Ain b upast m

load sys.mat

Ts = 50e-3; % 50 ms sampling time
n_delay = [1;21]; % delay as multiple of sampling time
% n_delay = 0;

xsize = length(A);
usize = size(B,2);
ysize = size(C,1);
dsize = 2; % number of disturbances



p = 200; % prediction horizon
m = 8; % control horizon



%% Define augmented system
% delay of n_delay time steps in 2nd component of u
% constant output disturbance

Adelay1 = [zeros(n_delay(1)-1,1), eye(n_delay(1)-1); zeros(1,n_delay(1))]; % delay component of A
Adelay2 = [zeros(n_delay(2)-1,1), eye(n_delay(2)-1); zeros(1,n_delay(2))]; % delay component of A
Adist = eye(dsize); % disturbance component of A

if (n_delay(1)==0)
    Aaug = [A, B(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(2),xsize), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(2)), Adist];

    Baug = [B(:,1), zeros(xsize,1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1];
        zeros(dsize,2)];
else
    Aaug = [A, B(:,1), zeros(xsize,n_delay(1)-1), B(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(1),xsize), Adelay1, zeros(n_delay(1),n_delay(2)), zeros(n_delay(1),dsize);
        zeros(n_delay(2),xsize), zeros(n_delay(2),n_delay(1)), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(1)), zeros(dsize,n_delay(2)), Adist];

    Baug = [zeros(xsize,usize);
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
orig_xsize = xsize;
xsize = length(Aaug);

%% Design observer
% expectation of output disturbance
% Qn = eye(orig_xsize+dsize);
% Qn = [eye(orig_xsize), zeros(orig_xsize,dsize);
%     zeros(dsize,orig_xsize), 1*eye(dsize)];
% Rn = zeros(dsize);

Qn = eye(orig_xsize);
Rn = eye(dsize);

% state disturbance to state matrix
% G = [eye(orig_xsize),zeros(orig_xsize,dsize);
%     zeros(n_delay,orig_xsize+dsize);
%     zeros(dsize,orig_xsize), eye(dsize)];

G = [zeros(xsize-dsize,orig_xsize);
    zeros(dsize,1), eye(dsize)*10, zeros(dsize,orig_xsize-dsize-1)];

Hkalman = [zeros(dsize,orig_xsize-dsize), eye(dsize)];

sys_kalman = ss(A,[B G], C, [D Hkalman], Ts);

[KEST,L,P,M,Z] = kalman(sys_kalman,Qn,Rn); % get observer values

%% Define QP matrices

YW=diag([1e1 1e3]');
UW=diag([50 1]');
% UW = diag([50 100]');
% YW = diag([1e4 1e5]');

YWT = kron(eye(p),YW);
UWT = kron(eye(m),UW);

% Y = Su*U + Sx*X
Su = zeros(ysize*p,usize*m);

for i=1:p
    for j=1:i
        % for first m inputs, make new columns
        if j<=m
            Su(1+(i-1)*ysize:i*ysize,1+(j-1)*usize:j*usize) = C*A^(i-j)*B;
            
        % m+1:p inputs are the same as input m
        else
            Su(1+(i-1)*ysize:i*ysize,1+(m-1)*usize:m*usize) = Su(1+(i-1)*ysize:i*ysize,1+(m-1)*usize:m*usize) +  C*A^(i-j)*B;
        end
    end
end


Sx = zeros(ysize*p,xsize);

for i=1:p
    Sx(1+(i-1)*ysize:i*ysize,:) = C*A^i;
end

H = Su'*YWT*Su + UWT;

Ga = YWT*Su;
Gb = Sx'*YWT*Su;
Gc = Su'*YWT*Su;


%% Define upper/lower bounds

% lb = [-0.3; 0];
% ub = [0.3; 1];
lb = [-0.1; 0];
ub = [0.1; 1];
LB = repmat(lb,m,1);
UB = repmat(ub,m,1);

% lb_rate = [0.1; 0.1];
% ub_rate = [0.1; 1];
lb_rate = [0.05;0.001];
ub_rate = [0.05;0.5];
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
  
xinit = zeros(xsize,1); 
upast = zeros(usize,1);

disp('MPC Problem Formulated!');

generate_file_v3m;

disp('Embedded files generated.')

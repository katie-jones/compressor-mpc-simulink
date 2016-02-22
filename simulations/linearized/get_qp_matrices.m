function [A,B,C,H,Ga,Gb,Gc] = get_qp_matrices(xinit,upast)
%#eml
%% Constants
% p = 200;
% m = 8;
p = 100;
m = 2;
Ts = 0.05;

%% Linearized system
[Ac,Bc,Cc] = get_linearized_matrices(xinit,upast);
[Ainit,Binit,Cinit] = discretize_rk4(Ac,Bc,Cc);
% Ainit = [
%     0.9351    0.0189   -0.0440   -0.0000    0.0449
%     0.0466    0.9200    0.0866    0.0000   -0.0884
%     0.9348   -0.7471    0.8074    0.0009    0.0598
%    -4.7533    3.8001   -8.8848    0.9907   -0.2001
%          0         0         0         0    0.9048];
%      
% Binit = [
%        -0.0001    0.0017
%     0.0001   -0.0032
%     0.0032    0.0015
%     7.0008   -0.0036
%          0    0.0675
%          ];
%      
% Cinit = [
%              0    1.0000         0         0         0
%    25.1590  -20.0646  100.0000         0         0
%    ];

n_delay = [0;40]; % delay as multiple of sampling time
% n_delay = 0;
% 
% xsize = length(Ainit);
% usize = size(Binit,2);
% ysize = size(Cinit,1);

xsize = 5;
usize = 2;
ysize = 2;
dsize = 2; % number of disturbances

%% Define augmented system
% delay of n_delay time steps in 2nd component of u
% constant output disturbance

Adelay2 = [zeros(n_delay(2)-1,1), eye(n_delay(2)-1); zeros(1,n_delay(2))]; % delay component of A
Adist = eye(dsize); % disturbance component of A

if (n_delay(1)==0)
    A = [Ainit, Binit(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(2),xsize), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(2)), Adist];

    B = [Binit(:,1), zeros(xsize,1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1];
        zeros(dsize,2)];
else
    Adelay1 = [zeros(n_delay(1)-1,1), eye(n_delay(1)-1); zeros(1,n_delay(1))]; % delay component of A

    A = [Ainit, Binit(:,1), zeros(xsize,n_delay(1)-1), Binit(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(1),xsize), Adelay1, zeros(n_delay(1),n_delay(2)), zeros(n_delay(1),dsize);
        zeros(n_delay(2),xsize), zeros(n_delay(2),n_delay(1)), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(1)), zeros(dsize,n_delay(2)), Adist];

    B = [zeros(xsize,usize);
        [zeros(n_delay(1)-1,1); 1], zeros(n_delay(1),1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1); 1];
        zeros(dsize,usize)];
end

% disturbance to output matrix
Cd = eye(ysize,dsize); 

C = [Cinit, zeros(ysize,sum(n_delay)), Cd];

%% Define correct global variables
xsize = length(A);
% xsize = 29;

%% Define weight matrices

YW=diag([1e1 1e3]');
UW=diag([50 1]');
% UW = diag([50 100]');
% YW = diag([1e4 1e5]');

YWT = kron(eye(p),YW);
UWT = kron(eye(m),UW);

%% Define system matrices

% Y = Su*U + Sx*X
Su = zeros(ysize*p,usize*m);

for i=1:p
    for j=1:i
        % for first m inputs, make new columns
        if j<=m
            Su(1+(i-1)*ysize:i*ysize,1+(j-1)*usize:j*usize) = C*A^(i-j)*B;
            
        % m+1:p inputs are the same as input m
        else
            toadd = C*A^(i-j)*B;
            for k=1:ysize
                Su(k+(i-1)*ysize,1+(m-1)*usize:m*usize) = Su(k+(i-1)*ysize,1+(m-1)*usize:m*usize) + toadd(k,:);
            end
%             Su(1+(i-1)*ysize:i*ysize,1+(m-1)*usize:m*usize) = Su(1+(i-1)*ysize:i*ysize,1+(m-1)*usize:m*usize) +  C*A^(i-j)*B;
%             Su(1+(i-1)*2:i*2,1+(8-1)*2:8*2) = Su(1+(i-1)*2:i*2,1+(8-1)*2:8*2) +  C(1:2,1:29)*A(1:29,1:29)^(i-j)*B(1:29,1:2);
        end
    end
end

% for i=m+1:p
%     offset = (i-1)*usize;
%     Su(:,1+(m-1)*usize:m*usize) = Su(:,1+(m-1)*usize:m*usize) + Su(:,1+offset:offset+usize);
% end


Sx = zeros(ysize*p,xsize);
Sx(1:ysize,:) = C*A;
for i=2:p
    toadd2 = C*A^i;
    for j=1:ysize
        Sx(j+(i-1)*ysize,:) = Sx(j+(i-2)*ysize,:) + toadd2(j,:);
    end
%     Sx(1+(i-1)*ysize:i*ysize,:) = Sx(1+(i-2)*ysize:(i-1)*ysize,:) + C*A^i;
end

%% Calculate QP matrices

H = Su'*YWT*Su + UWT;

Ga = YWT*Su;
Gb = Sx'*YWT*Su;
Gc = Su'*YWT*Su;

end
clear all; close all; %clc
% cd qpOASES-3.0beta\interfaces\matlab\ 
% path = [path pwd];
% cd ..\..\..
addpath('qpOASES-3.0beta\interfaces\matlab\')
global xinit H A B C M ysize xsize usize p LB UB Ga Gb Gc Ain b m upast

load sys
Ts = 0.05;
% delay free system
% [A,B,C,D,Ts]=ssdata(sys);
sys = ss(A,B,C,D,Ts);


[xsize dmmy]=size(A);
[dmmy usize]=size(B);
[ysize dmmy]=size(C);

MaxDelay = 40;% max delay as a multiple of sampling time
MaxDelay = 3;

Isize = usize*MaxDelay;

Adelay = [zeros(usize,Isize+usize);eye(Isize) zeros(Isize,usize)];

Bdelay = [eye(usize,usize);zeros(Isize,usize)];

IDM = [0;20]; % input to state delay array (or matrix) as a multiple of sampling time 
IDM = [0;2];
Edelay = zeros(xsize,Isize+usize);

% Additional for loop needed to individualy insert input to state delay
for i=1:usize
        Edelay(:,IDM(i,1)*usize+i)=B(:,i);
end

Aaug = [A Edelay;zeros(Isize+usize,xsize) Adelay];
Baug = [zeros(xsize,usize);Bdelay];
Caug = [C zeros(ysize,Isize+usize)];

sysa = ss(Aaug,Baug,Caug,D,Ts);




[A,B,C,D]=ssdata(sysa);



[xsize dmmy]=size(A);
[dmmy usize]=size(B);
[ysize dmmy]=size(C);

%Bd=B(:,1);
Bd=zeros(xsize,1);


Dd=ones(ysize,1);


Abar=[1 0;0 1];
Bbar=[10 0;0 10];
Cbar=[1 1];
Dbar=[1 1];

[xdsize dmmy]=size(Abar);


A1=A;
A2=Bd*Cbar*0;
A4=Abar;

A3=zeros(xdsize,xsize);

A=[[A1 A2];[A3 A4]];


B=[B;zeros(xdsize,usize)];
BD=[Bd*Dbar;Bbar];
BN=zeros(xsize+xdsize,ysize);

BDN=[BD BN];

BO=0*B(:,1);

%C=[C Dd*Cbar];
C=[C eye(2)];

DO=0*D(:,1);
Dbarm=zeros(ysize,2);
Dnoise=eye(ysize);

sys=ss(A,[B BO BDN],C,[D DO Dbarm Dnoise],10);

[KEST,L,W,M,Z] = kalman(sys,eye(5),eye(2));

global xinit H A B C M ysize xsize usize p LB UB Ga Gb Gc Ain b m upast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [A,B,C,D] system model
% P         the prediction model matrix
% H         weight matrix for QP u'*H*u + f'*u
% xsize     system state size
% ysize     system output size
% usize     system input size
% p         prediction horizon
% m         move horizon
% LB        lower bounds for QP (only for inputs)
% UB        upper bounds for QP (only for inputs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtain system sizes

xsize=size(A,1);
usize=size(B,2);
ysize=size(C,1);

% prediction and move horizons

p=200;
m=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prediction matrix P

P=zeros(p*ysize,p*usize+xsize);

for i=1:p
    P(1+(i-1)*ysize:ysize+(i-1)*ysize,1:xsize)=C*A^i;
end


for j=1:p
    for i=1:p-j+1
        P(1+(i-1)*ysize+(j-1)*ysize:(j-1)*ysize+ysize+(i-1)*ysize,...
        xsize+1+(j-1)*usize:xsize+usize*j)=C*A^(i-1)*B;
    end
end


for i=m:p-1
    P(i*ysize+1:(i+1)*ysize,xsize+(m-1)*usize+1:xsize+m*usize)=...
    P(i*ysize+1:(i+1)*ysize,xsize+(m-1)*usize+1:xsize+m*usize)+...
    P((i-1)*ysize+1:i*ysize,xsize+(m-1)*usize+1:xsize+m*usize);
end
   
P=P(:,1:xsize+m*usize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output weights for QP

YW=diag([1e1 1e3]');

for i=1:p
    YWT(1+(i-1)*ysize:i*ysize,1+(i-1)*ysize:i*ysize)=YW;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input weights for QP

UW=diag([50 1]'); % UW=diag([50 5]');

UWT=zeros(m*usize);

for i=1:m
    UWT(1+(i-1)*usize:i*usize,1+(i-1)*usize:i*usize)=UW;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upper bounds for the inputs

% UBU=[0 0]'; 
UBU=[0.3 1]';

UB=inf*ones(m*usize,1);

for i=1:m
    UB(1+(i-1)*usize:i*usize)=UBU;
end

% lower bounds for the inputs

% LBU=[0 0]'; 
LBU=[-0.3 0]';

LB=-inf*ones(m*usize,1);

for i=1:m
    LB(1+(i-1)*usize:i*usize)=LBU;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ueye=eye(usize);

% input rate constraints prev labels refer to previous time step 

in_ratedown = [0.05;0.001]; % in_ratedown = [0.05;0.5];

in_rateup   = [0.05;0.5]; % in_rateup   = [0.05;1];

for i=1:m
    b_down((i-1)*usize+1:i*usize,1)=in_ratedown;
end

for i=1:m
    b_up((i-1)*usize+1:i*usize,1)=in_rateup;
end

b=[b_down;b_up];


A_RATEDOWN=zeros((m-1)*usize,m*usize);

for i=1:m-1
    A_RATEDOWN((i-1)*usize+1:i*usize,(i-1)*usize+1:i*usize+usize)=[ueye -ueye];
end

A_DOWNprev=zeros(usize,m*usize);
A_DOWNprev(1:usize,1:usize)=-eye(usize);


A_RATEUP=zeros((m-1)*usize,m*usize);

for i=1:m-1
    A_RATEUP((i-1)*usize+1:i*usize,(i-1)*usize+1:i*usize+usize)=[-ueye ueye];
end

A_UPprev=zeros(usize,m*usize);
A_UPprev(1:usize,1:usize)=eye(usize);


Ain=[A_DOWNprev;A_RATEDOWN;A_UPprev;A_RATEUP];


% Hessian matrix and Gradient for the QP

P1=P(:,1:xsize);
P2=P(:,xsize+1:m*usize+xsize);

H=P2'*YWT*P2+UWT;

Ga=YWT*P2;
Gb=P1'*YWT*P2;
Gc=P2'*YWT*P2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize system state

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

% Function to linearize discharge tank about point xbar, ubar

function [A, B, C] = linearize_tank(x,u)
%#eml

%% Constants

xsize = 5;
usize = 4;
ysize = 2;

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);
pd = x(2*xsize+1);

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

VolumeT = const_tank();
[~,~,~,C_coeff,~,~,D2] = comp_coeffs();

[SpeedSound,Pin,Pout,V1,V2] = const_flow();

ud1 = u1(3);
uin2 = u2(2);

pin_tank = x1(2);
pout_tank = x2(1);



%% Dynamics of discharge tank
Att = - getValveDerivative(pd,pout_tank,SpeedSound,VolumeT,uin2,C_coeff) - getValveDerivative(pin_tank,pd,SpeedSound,VolumeT,ud1,D2);

%% Interaction between compressors and discharge tanks

% Effect of compressors on tank
Act = [0; getValveDerivative(pin_tank,pd,SpeedSound,VolumeT,ud1,D2); 0; 0; 0;
    getValveDerivative(pd,pout_tank,SpeedSound,VolumeT,uin2,C_coeff); 0; 0; 0; 0]';

% Effect of tank on compressors
Atc = [0; getValveDerivative(pin_tank,pd,SpeedSound,V2,ud1,D2); 0; 0; 0;
    getValveDerivative(pd,pout_tank,SpeedSound,V1,uin2,C_coeff); 0; 0; 0; 0];

%% Compressor dynamics
[A1,B1,C1] = get_linearized_matrices(x1,[u1; Pin; pd]);

[A2,B2,C2] = get_linearized_matrices(x2,[u2; pd; Pout]);


%% Output system
A = [ [A1, zeros(xsize);
    zeros(xsize), A2;
    Act], [Atc; Att] ];

B = [B1, zeros(xsize,2);
    zeros(xsize,2), B2;
    zeros(1,2*2)];

% Outputs: p21,SD1,p22,SD2,PD
% C = [C1, zeros(ysize,xsize), zeros(ysize,1);
%     zeros(ysize,xsize), C2, zeros(ysize,1);
%     zeros(1,2*xsize), 1];

% Outputs: p21,SD1,p22,SD2,PD
C = [C1, zeros(ysize,xsize), zeros(ysize,1);
    zeros(ysize,xsize), C2, zeros(ysize,1);
    [0,0,1,0,0],[0,0,-1,0,0],0];



end



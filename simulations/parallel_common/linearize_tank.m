% Function to linearize discharge tank about point xbar, ubar

function [A, B, C] = linearize_tank(x,u)
%#eml

%% Constants
pd = x(end);
ud = u(end);

xsize = 5;
usize = 6;
ysize = 2;

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

u1(end) = pd;
u2(end) = pd;

[~,VolumeT,D2_t] = const_tank();
[~,~,~,~,~,~,D2] = comp_coeffs();

[SpeedSound,~,~,~,V2] = const_flow();

ud1 = u1(3);

pin_tank = x1(2);
pout_tank = x2(1);



%% Dynamics of discharge tank
Att = - getValveDerivative(pd,pout_tank,SpeedSound,VolumeT,ud,D2_t) - getValveDerivative(pin_tank,pd,SpeedSound,VolumeT,ud1,D2);

%% Interaction between compressors and discharge tanks

% Effect of compressors on tank
Act = [0; getValveDerivative(pin_tank,pd,SpeedSound,VolumeT,ud1,D2); 0; 0; 0;
    getValveDerivative(pd,pout_tank,SpeedSound,VolumeT,ud,D2); 0; 0; 0; 0]';

% Effect of tank on compressors
Atc = [0; getValveDerivative(pin_tank,pd,SpeedSound,V2,ud1,D2); 0; 0; 0;
    0; getValveDerivative(pd,pout_tank,SpeedSound,V2,ud,D2); 0; 0; 0];

%% Compressor dynamics
[A1,B1,C1] = get_linearized_matrices(x1,u1);

[A2,B2,C2] = get_linearized_matrices(x2,u2);


%% Output system
A = [ [A1, zeros(xsize);
    zeros(xsize), A2;
    Act], [Atc; Att] ];

B = [B1, zeros(xsize,2);
    zeros(xsize,2), B2;
    zeros(1,2*2)];

% Outputs: p21,SD1,p22,SD2,PD
C = [C1, zeros(ysize,xsize), zeros(ysize,1);
    zeros(ysize,xsize), C2, zeros(ysize,1);
    zeros(1,2*xsize), 1];



end


function dp2 = getValveDerivative(p1,p2,SpeedSound,V,u,D)
dp2 = SpeedSound*SpeedSound/V * 1e-5 * 100/2/sqrt(abs(p1*100-p2*100)) * [u^3, u^2, u, 1] * D(1:4)';

end


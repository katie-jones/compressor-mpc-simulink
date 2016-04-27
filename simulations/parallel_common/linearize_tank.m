% Function to linearize discharge tank about point xbar, ubar

function [A, B, C] = linearize_tank(x,u)
%#eml

%% Constants

xsize = 5;
usize = 6;
ysize = 2;

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

[~,~,~,~,~,~,D2] = comp_coeffs();

[SpeedSound,~,~,Vin,Vout] = const_flow();

ud1 = u1(3);

%% Interaction between compressors

dc1c2 = get_Act(x1(2),ud1,x2(1),SpeedSound,Vin,D2); % effect of comp 1 on comp 2
dc2c1 = get_Act(x1(2),ud1,x2(1),SpeedSound,Vout,D2); % effect of comp 2 on comp 1

Ac1c2 = zeros(xsize);
Ac2c1 = Ac1c2;

Ac1c2(2,1) = dc1c2;
Ac2c1(1,2) = dc2c1;


%% Compressor dynamics
[A1,B1,C1] = get_linearized_matrices(x1,u1);

[A2,B2,C2] = get_linearized_matrices(x2,u2);


%% Output system
A = [ A1, Ac2c1; Ac1c2, A2 ];

B = [B1, zeros(xsize,2);
    zeros(xsize,2), B2];


% Outputs: p21,p22,SD1,SD2
C = [C1, zeros(ysize,xsize);
    zeros(ysize,xsize), C2];



end

% Get component of A matrix that is the effect of compressor pressure on
% adjoining tank
function Act = get_Act(p2,ud,pd,SpeedSound,VolumeT,D2)

Act = SpeedSound*SpeedSound/VolumeT * 1e-5 * (1/2*100/sqrt(abs(p2*100-pd*100))) * (D2(1)*ud^3 + D2(2)*ud^2 + D2(3)*ud + D2(4));
    
end


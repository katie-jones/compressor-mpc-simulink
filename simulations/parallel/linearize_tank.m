% Function to linearize discharge tank about point xbar, ubar

function [A, B, C] = linearize_tank(x,u)

%% Constants
pd = x(end);
udt = u(end);

xsize = 5;
usize = 5;
ysize = 2;

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

[Out_pres_t,VolumeT,D2_t] = tank_params();
[~,~,~,~,~,~,D2] = comp_coeffs();

[SpeedSound,~,~,V1,V2] = flow_params();

ud1 = u1(3);
ud2 = u2(3);




%% Dynamics of discharge tank
Att = -SpeedSound * SpeedSound / VolumeT * 1e-5 * (1/2*100/sqrt(abs(pd*100-Out_pres_t*100))) * (D2_t(1)*udt^3 + D2_t(2)*udt^2 + D2_t(3)*udt + D2_t(4));



%% Interaction between compressors and discharge tanks

Ac1t = get_Act(x1,udt,pd,SpeedSound,VolumeT,D2_t);

Ac2t = get_Act(x2,udt,pd,SpeedSound,VolumeT,D2_t);

Atc1 = get_Atc(x1,ud1,pd,SpeedSound,V1,D2);

Atc2 = get_Atc(x2,ud2,pd,SpeedSound,V2,D2);

%% Compressor dynamics
[A1,B1,C1] = get_linearized_matrices(x1,u1);

[A2,B2,C2] = get_linearized_matrices(x2,u2);


%% Output system
A = [A1, zeros(xsize);
    zeros(xsize), A2;
    Ac1t, Ac2t];

A = [A, [Atc1; Atc2; Att]];

B = [B1, zeros(xsize,2);
    zeros(xsize,2), B2;
    zeros(1,2*2)];

% Outputs: SD1, SD2, p21-p22, pd
% C = [C1(2,:), zeros(1,xsize), 0;
%     zeros(1,xsize), C2(2,:), 0;
%     C1(1,:), -C2(1,:), 0;
%     zeros(1,2*xsize), 1];
C = [C1, zeros(ysize,xsize), zeros(ysize,1);
    zeros(ysize,xsize), C2, zeros(ysize,1);
    zeros(1,2*xsize), 1];



end

% Get component of A matrix that is the effect of compressor states on pd
function Act = get_Act(x,ud,pd,SpeedSound,VolumeT,D2)

p2 = x(2);

Act = [0, ...
    SpeedSound*SpeedSound/VolumeT * 1e-5 * (1/2*100/sqrt(abs(p2*100-pd*100))) * (D2(1)*ud^3 + D2(2)*ud^2 + D2(3)*ud + D2(4)),...
    0, 0, 0
    ];


end


% Get component of A matrix that is the effect of pd on compressor states
function Atc = get_Atc(x,ud,pd,SpeedSound,VolumeT2,D2)

p2 = x(2);

Atc = [0;
    SpeedSound * SpeedSound / VolumeT2 * 1e-5 * (1/2*100/sqrt(abs(p2*100-pd*100))) * (D2(1)*ud^3 + D2(2)*ud^2 + D2(3)*ud + D2(4));
    0; 0; 0];

end


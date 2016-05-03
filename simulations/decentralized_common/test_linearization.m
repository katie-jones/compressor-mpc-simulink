function [systf,rise_time] = test_linearization()

volumes=2.5:5:25;
n = length(volumes);
systf = cell(n,1);
rise_time = zeros(n,1);

x_init_lin = [0.916 1.145 0.152 440 0]';
pd = 1.12;

uoff = [0.304 0.43 1 0 0]';
ud = 0.7;

x = [x_init_lin; x_init_lin; pd];
u = [uoff; uoff; ud];
legstr = cell(n,1);

figure; hold on
for i=1:n
    [A,B,C] = linearize_tank(x,u,volumes(i));
    systf{i} = tf(ss(A,B,C,0));
    [y,t] = step(systf{i}(3,1));
    S = stepinfo(y,t);
    rise_time(i) = S.RiseTime;
    step(systf{i}(4,1));%,systf{i}(1,1))
    legstr{i} = ['Volume: ',num2str(volumes(i))];
end
legend(legstr,'location','northeast')
title('Step response: torque input to opposite compressor surge distance')
ylabel('Surge distance')
figure; hold on
for i=1:n
    bode(systf{i}(4,1))
end
legend(legstr,'location','southwest')

end


% Function to linearize discharge tank about point xbar, ubar

function [A, B, C] = linearize_tank(x,u,VolumeT)


%% Constants
pd = x(end);
ud = u(end);

xsize = 5;
usize = 5;
ysize = 2;

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

u1(end) = pd;
u2(end) = pd;

[Out_pres_t,~,D2_t] = const_tank();
[~,~,~,~,~,~,D2] = comp_coeffs();

[SpeedSound,~,~,~,V2] = const_flow();

ud1 = u1(3);
ud2 = u2(3);




%% Dynamics of discharge tank
Att = -getValveDerivative(pd,Out_pres_t,SpeedSound,VolumeT,ud,D2_t) - getValveDerivative(x1(2),pd,SpeedSound,VolumeT,ud1,D2) - getValveDerivative(x2(2),pd,SpeedSound,VolumeT,ud2,D2);

%% Interaction between compressors and discharge tanks

% Effect of compressors on tank
Act = [0; getValveDerivative(x1(2),pd,SpeedSound,VolumeT,ud1,D2); 0; 0; 0;
    0; getValveDerivative(x2(2),pd,SpeedSound,VolumeT,ud2,D2); 0; 0; 0]';

% Effect of tank on compressors
Atc = [0; getValveDerivative(x1(2),pd,SpeedSound,V2,ud1,D2); 0; 0; 0;
    0; getValveDerivative(x2(2),pd,SpeedSound,V2,ud2,D2); 0; 0; 0];

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


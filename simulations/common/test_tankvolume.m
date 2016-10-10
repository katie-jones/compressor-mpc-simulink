function test_tankvolume()

Vtank = 0.5:2:8.5;
x_init_lin = [0.916; 1.145; 0.152; 440; 0];
P_D = 1.12;
[Ts, xsize_comp, xsize, usize_comp, ysize, uoff1, uoff2, ud] = const_sim();

xinit = [x_init_lin; x_init_lin; P_D];
uinit = [uoff1; uoff2; ud];

t_plot = linspace(0,100,100);

y = zeros(length(t_plot),length(Vtank));
ycross = y;

legstr = cell(length(Vtank),1);

for i=1:length(Vtank)
    [A,B,C] = linearize_tank(xinit, uinit, Vtank(i));
    sys{i} = ss(A,B,C,0);
    [y_in,t_in] = step(sys{i});
    y(:,i) = interp1(t_in,squeeze(y_in(:,2,1)),t_plot);
    ycross(:,i) = interp1(t_in,squeeze(y_in(:,2,3)),t_plot);
    legstr{i} = num2str(Vtank(i));
end

figure; 
hold on
plot(t_plot,y)
plot(t_plot,-ycross,'linestyle','--')
legend(legstr)
title('Step response to torque input')
xlabel('Time [s]')
ylabel('Surge Distance')

xlim([0 50])



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


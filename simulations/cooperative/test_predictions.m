% Demonstrate effect of Td on surge distances for various conditions

addpath('../common')
addpath('../parallel_common/')

save_plots=1;

%% Constants
[Ts, xsize_comp, xsize, ~, ysize] = const_sim();
[n_delay,dsize,ucontrolsize,p,m] = const_mpc();

xtotalsize = xsize + 2*sum(n_delay) + 2*dsize;

P_D = 1.12;

x_init_lin = [0.916; 1.145; 0.152; 440; 0];

xinit = [x_init_lin; x_init_lin; P_D];

uinit = zeros(2*ucontrolsize,1);

%% Predictions
[A,B,C,H1,H2,Ga1,Ga2,Gb1,Gb2,Gc1,Gc2,dx,Sx,Gd1,Gd2,Sf,Su1,Su2] = get_qp_matrices(xinit,uinit);

du = repmat([0.01; 0],m,1); % step on torque, no recycle valve

y1 = reshape(Su1*du,ysize,p)'; % prediction of compressor 1

t = 1:p;

SD_sum = zeros(p,2);

SD_sum(1,:) = [0 0];
for i=2:p
    SD_sum(i,:) = SD_sum(i-1,:) + Ts*y1(i,1:2);
end
%%

fig=figure; 
set(fig,'units','centimeters','position',[35 5 28 15])
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',2)
subplot(1,2,1)
plot(t,y1(:,1:2));
grid on
title('Predicted surge distance')
xlabel({'Horizon length (m)'; 'x = [0.92 1.15 0.15 440 0]'''})
ylabel('Relative surge distance [%]')
legend('SD_1','SD_2','location','southwest')

subplot(1,2,2)
plot(t,SD_sum.^2);

ylim([0 0.25])
xlim([0 p])
grid on
title('Cumulative effect of predicted surge distance')
xlabel({'Horizon length (m)'; '\Delta T_{d,1}=0.01, p=2'})
ylabel('Cumulative effect of surge distance')

if save_plots
    fig=printplot(fig);
    saveas(fig,'Td_effect.pdf')
end

set(0,'defaultlinelinewidth',1)


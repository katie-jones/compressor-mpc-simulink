%% Plot settings
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',1)

fig=figure;
% set(fig,'units','centimeters','position',[10 10 10 7])
plot(pout.time,pout.signals.values+1.126);
title('Compressor outlet pressure')
xlabel('Time [s]')
ylabel('Pressure [bar?]')
leg=legend('Compressor 1','Compressor 2');
set(leg,'fontsize',14,'location','best')


fig=figure;
plot(PD.time,PD.signals.values); hold on
plot([min(PD.time),max(PD.time)],[yref(4),yref(4)],':k')
title('Tank pressure')
xlabel('Time [s]')
ylabel('Pressure [bar?]')
legend('Value','Reference')

fig=figure; 
plot(SD.time,SD.signals.values);
hold on
plot([min(SD.time),max(SD.time)],[yref(1),yref(1)],':k')
title('Surge distance')
xlabel('Time [s]')
ylabel('Relative surge distance [%]')
leg=legend('Compressor 1','Compressor 2','Reference');
set(leg,'fontsize',14,'location','southeast')

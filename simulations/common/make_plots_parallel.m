%% Plot settings
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

fig=figure;
% set(fig,'units','centimeters','position',[10 10 10 7])
plot(pout.time,pout.signals.values+1.126);
grid on
title('Compressor outlet pressure')
xlabel('Time [s]')
ylabel('Pressure [bar?]')
leg=legend('Compressor 1','Compressor 2');
set(leg,'fontsize',14,'location','best')

if save_plots
    fig=printplot(fig);
    saveas(fig,[resultsdir,'/pout.pdf'])
end


fig=figure;
plot(PD.time,PD.signals.values); hold on
plot([min(PD.time),max(PD.time)],[yref(4),yref(4)],':k')
grid on
title('Tank pressure')
xlabel('Time [s]')
ylabel('Pressure [bar?]')
leg = legend('Value','Reference');
set(leg,'fontsize',14,'location','best')
if save_plots
    fig=printplot(fig);
    saveas(fig,[resultsdir,'/pd.pdf'])
end


fig=figure; 
plot(SD.time,SD.signals.values);
hold on
grid on
plot([min(SD.time),max(SD.time)],[yref(1),yref(1)],':k')
title('Surge distance')
xlabel('Time [s]')
ylabel('Relative surge distance [%]')
leg=legend('Compressor 1','Compressor 2','Reference');
set(leg,'fontsize',14,'location','southeast')
if save_plots
    fig=printplot(fig);
    saveas(fig,[resultsdir,'/sd.pdf'])
end


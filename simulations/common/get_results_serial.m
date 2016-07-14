filename = '/home/katie/school/MasterThesis/cpp/serial/output/cent_output.dat';
fname_info = '/home/katie/school/MasterThesis/cpp/serial/output/cent_info.dat';
res_cent = read_cpp_results(filename,fname_info,10,4,4);

n_solver_iterations = 3;

filename = ['/home/katie/school/MasterThesis/cpp/serial/output/coop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/serial/output/coop_info',num2str(n_solver_iterations),'.dat'];
res_coop = read_cpp_results(filename,fname_info,10,4,4);

filename = ['/home/katie/school/MasterThesis/cpp/serial/output/ncoop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/serial/output/ncoop_info',num2str(n_solver_iterations),'.dat'];
res_ncoop = read_cpp_results(filename,fname_info,10,4,4);

td_offset = 0.304;

u_dist = [0.43, 0.43, 0.33, 0.33];
t_dist = [0, 50, 50, res_cent.t(end)];

%% Time responses
res_cent.sd = [res_cent.y(:,2) - res_cent.yref(2), res_cent.y(:,4) - res_cent.yref(4)];
res_coop.sd = [res_coop.y(:,2) - res_coop.yref(2), res_coop.y(:,4) - res_coop.yref(4)];
res_ncoop.sd = [res_ncoop.y(:,2) - res_ncoop.yref(2), res_ncoop.y(:,4) - res_ncoop.yref(4)];

res_cent.p = res_cent.y(:,[1,3]);
res_coop.p = res_coop.y(:,[1,3]);
res_ncoop.p = res_ncoop.y(:,[1,3]);

res_cent.td = res_cent.u(:,[1,3]) + td_offset;
res_coop.td = res_coop.u(:,[1,3]) + td_offset;
res_ncoop.td = res_ncoop.u(:,[1,3]) + td_offset;

res_cent.ur = res_cent.u(:,[2,4]);
res_coop.ur = res_coop.u(:,[2,4]);
res_ncoop.ur = res_ncoop.u(:,[2,4]);

uwt = 1e2*[1 10 1 10]';
ywt = [1 10 1 10]';
res_cent = AddCostFunctionSerial(res_cent,uwt,ywt);
res_coop = AddCostFunctionSerial(res_coop,uwt,ywt);
res_ncoop = AddCostFunctionSerial(res_ncoop,uwt,ywt);

results = {res_cent, res_coop, res_ncoop};

figs1 = PlotResults(results,'p');

xlims = [0, 320];

for i=1:length(figs1)
    fig1=figs1{i};
    figure(fig1);
    title('Tank Output Pressure')
    xlabel('Time [s]')
    ylabel('Pressure [atm]')
    legend('Centralized','Cooperative','Non-cooperative')
    xlim(xlims);
    
    if exist('saveplots','var') && (saveplots ~= 0)
        saveas(fig1,['serial_p',num2str(i),'.fig']);
        saveas(fig1,['serial_p',num2str(i),'.pdf']);
    end
end


figs1 = PlotResults(results,'sd');
for i=1:length(figs1)
    fig1 = figs1{i};
    figure(fig1);
    title('Surge Distance')
    xlabel('Time [s]')
    ylabel('Relative Surge Control Distance [%]')
    xlim(xlims);
    
    if exist('saveplots','var') && (saveplots ~= 0)
        saveas(fig1,['serial_sd',num2str(i),'.fig']);
        saveas(fig1,['serial_sd',num2str(i),'.pdf']);
    end
end


figs1 = PlotResults(results,'td');
for i=1:length(figs1)
    fig1 = figs1{i};
    figure(fig1);
    title('Normalized Torque Setting')
    xlabel('Time [s]')
    ylabel('Torque')
    xlim(xlims);
    
    if exist('saveplots','var') && (saveplots ~= 0)
        saveas(fig1,['serial_td',num2str(i),'.fig']);
        saveas(fig1,['serial_td',num2str(i),'.pdf']);
    end
end

figs1 = PlotResults(results,'ur');
for i=1:length(figs1)
    fig1 = figs1{i};
    figure(fig1);
    title('Normalized Recycle Valve Opening')
    xlabel('Time [s]')
    ylabel('Valve Opening')
    xlim(xlims);
    
    if exist('saveplots','var') && (saveplots ~= 0)
        saveas(fig1,['serial_ur',num2str(i),'.fig']);
        saveas(fig1,['serial_ur',num2str(i),'.pdf']);
    end
end


% fig1 = figure;
% plot(t_dist,u_dist)
% title('Tank Output Valve Opening')
% ylabel('Valve Opening')
% xlabel('Time [s]')
% ylim([0, max(u_dist)*1.1])

%% Cost function

fig1 = PlotResults(results,'J');
title('Cost Function')
xlabel('Time [s]')
ylabel('Cost Function')
legend('Centralized','Cooperative','Non-cooperative')

if exist('saveplots','var') && (saveplots ~= 0)
    saveas(fig1{1},'serial_j.fig');
    saveas(fig1{1},'serial_j.pdf');
end

%% Cost function
% close all
% xlims = [40 140];
% figure; plot(res_ncoop.t,res_ncoop.y(:,[1,3]));
% hold on;
% plot(res_cent.t,res_cent.y(:,[1,3]),'--'); title('Pressures')
% xlim(xlims);
% figure; plot(res_ncoop.t,res_ncoop.y(:,[2,4]),res_cent.t,res_cent.y(:,[2,4])); title('Surge Distance')
% xlim(xlims);
% figure; plot(res_ncoop.t,res_ncoop.u(:,[2,4]),res_cent.t,res_cent.u(:,[2,4])); title('Recycle Valve Opening')
% xlim(xlims);
% figure; plot(res_ncoop.t,res_ncoop.u(:,[1,3]),res_cent.t,res_cent.u(:,[1,3])); title('Torque Input')
% xlim(xlims);

%%
% close all
%
% res_cent = AddCostFunctionSerial(res_cent);
% res_coop = AddCostFunctionSerial(res_coop);
% res_ncoop = AddCostFunctionSerial(res_ncoop);
%
% results = {res_cent, res_ncoop};
%
% inds = res_cent.t>=50.0;
% xlims = [50 300];
%
% titles = {'Td1','UR1','Td2','UR2';
%     'Po1','SD1','Po2','SD2'};
%
% fig=PlotCostFunctionSerial(results,'Ju',titles,0);

% fig_width=0.4;
% fig_height=0.3;
%
% for i=1:4
%     figure; plot(res_cent.t,[res_cent.Ju(:,i),res_coop.Ju(:,i)]);
%     figure; plot(res_cent.t,[res_cent.Jy(:,i),res_coop.Jy(:,i)]);
% end

% fig=PlotCostFunction(results,'Jsd','Surge Distance');
% set(fig,'units','normalized','position',[0,0.5,fig_height,fig_width])
% legend('Centralized','Cooperative (3 it.)','Non-cooperative (3 it.)')
% xlim(xlims);
% if saveplots
%     saveas(fig,'sd_penalty.fig');
%     saveas(fig,'sd_penalty.pdf');
% end
%
% fig=PlotCostFunction(results,'Jp','Output Pressure');
% set(fig,'units','normalized','position',[0.3,0.5,fig_height,fig_width])
% xlim(xlims);
% if saveplots
%     saveas(fig,'pout_penalty.fig');
%     saveas(fig,'pout_penalty.pdf');
% end
%
% fig=PlotCostFunction(results,'Jtd','Torque Setting');
% set(fig,'units','normalized','position',[0,0,fig_height,fig_width])
% xlim(xlims);
% if saveplots
%     saveas(fig,'td_penalty.fig');
%     saveas(fig,'td_penalty.pdf');
% end
%
% fig=PlotCostFunction(results,'Jur','Recycle Valve Opening');
% set(fig,'units','normalized','position',[0.3,0,fig_height,fig_width])
% xlim(xlims);
% if saveplots
%     saveas(fig,'ur_penalty.fig');
%     saveas(fig,'ur_penalty.pdf');
% end
% figure;
% plot(res_cent.t(inds)-50,[res_cent.Jsd(inds),res_coop.Jsd(inds),res_ncoop.Jsd(inds)]/norm_sd*100)
% title('Surge Distance')
% xlabel('Time [s]')
% ylabel('Penalty term')
% legend('Centralized','Cooperative','Non-Cooperative')
% xlim(xlims)
%
% figure;
% plot(res_cent.t(inds)-50,[res_cent.Jp(inds),res_coop.Jp(inds),res_ncoop.Jp(inds)]/norm_p*100)
% title('Output Pressure')
% xlabel('Time [s]')
% ylabel('Penalty term')
% xlim(xlims)
%
% figure;
% plot(res_cent.t(inds)-50,[res_cent.Jtd(inds),res_coop.Jtd(inds),res_ncoop.Jtd(inds)])
% title('Torque Setting')
% xlabel('Time [s]')
% ylabel('Penalty term')
%
% figure;
% plot(res_cent.t(inds)-50,[res_cent.Jur(inds),res_coop.Jur(inds),res_ncoop.Jur(inds)])
% title('Recycle Valve Opening')
% xlabel('Time [s]')
% ylabel('Penalty term')



% %%
% figure;
% plot(res_cent.t,[res_cent.y(:,1),res_coop.y(:,1),res_ncoop.y(:,1)]);
% legend('Centralized','Cooperative','Non-cooperative')
% title('Surge Distance Response')
%
% %%
% figure;
% plot(res_cent.t,[res_cent.y(:,4),res_coop.y(:,4),res_ncoop.y(:,4)]);
% legend('Centralized','Cooperative','Non-cooperative')
% title('Output Pressure Response')

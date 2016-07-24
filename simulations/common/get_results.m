filename = '/home/katie/school/MasterThesis/cpp/parallel/output/cent_output.dat';
fname_info = '/home/katie/school/MasterThesis/cpp/parallel/output/cent_info.dat';
res_cent = read_cpp_results(filename,fname_info,11,4,4);

n_solver_iterations = 3;

filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/coop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/parallel/output/coop_info',num2str(n_solver_iterations),'.dat'];
res_coop = read_cpp_results(filename,fname_info,11,4,4);

filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/noncoop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/parallel/output/noncoop_info',num2str(n_solver_iterations),'.dat'];
res_ncoop = read_cpp_results(filename,fname_info,11,4,4);

u_dist = [0.7, 0.7, 0.4, 0.4];
t_dist = [0, 50, 50, res_cent.t(end)];

td_offset = 0.304;

fname_base = 'parallel_';

%% Time responses

res_cent.sd = res_cent.y(:,1) - res_cent.yref(1);
res_coop.sd = res_coop.y(:,1) - res_coop.yref(1);
res_ncoop.sd = res_ncoop.y(:,1) - res_ncoop.yref(1);

res_cent.p = res_cent.y(:,4);
res_coop.p = res_coop.y(:,4);
res_ncoop.p = res_ncoop.y(:,4);

res_cent.td = res_cent.u(:,[1]) + td_offset;
res_coop.td = res_coop.u(:,[1]) + td_offset;
res_ncoop.td = res_ncoop.u(:,[1]) + td_offset;

res_cent.ur = res_cent.u(:,[2]);
res_coop.ur = res_coop.u(:,[2]);
res_ncoop.ur = res_ncoop.u(:,[2]);

uwt = 1e2*[1 1 10 10]';
ywt = [1 1 100]';
res_cent = AddCostFunction(res_cent,uwt,ywt);
res_coop = AddCostFunction(res_coop,uwt,ywt);
res_ncoop = AddCostFunction(res_ncoop,uwt,ywt);

results = {res_cent, res_coop, res_ncoop};

figs = PlotResults(results,'p');
title('Tank Output Pressure')
xlabel('Time [s]')
ylabel('Pressure [atm]')
legend('Centralized','Cooperative','Non-cooperative')

fname=[fname_base,'p'];
    
if exist('saveplots','var') && (saveplots ~= 0)
    saveas(figs{1},[fname,'.fig']);
    saveas(figs{1},[fname,'.pdf']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
    fname,'.tex'],'width','0.8\linewidth','figurehandle',figs{1},'showInfo',false);
end


figs = PlotResults(results,'sd');
title('Surge Distance')
xlabel('Time [s]')
ylabel('Relative Surge Control Distance [%]')

fname=[fname_base,'sd'];
    
if exist('saveplots','var') && (saveplots ~= 0)
    saveas(figs{1},[fname,'.fig']);
    saveas(figs{1},[fname,'.pdf']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
    fname,'.tex'],'width','0.8\linewidth','figurehandle',figs{1},'showInfo',false);
end


figs = PlotResults(results,'td');
title('Normalized Torque Setting')
xlabel('Time [s]')
ylabel('Torque')

fname=[fname_base,'td'];
    
if exist('saveplots','var') && (saveplots ~= 0)
    saveas(figs{1},[fname,'.fig']);
    saveas(figs{1},[fname,'.pdf']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
    fname,'.tex'],'width','0.8\linewidth','figurehandle',figs{1},'showInfo',false);
end


figs = PlotResults(results,'ur');
title('Normalized Recycle Valve Opening')
xlabel('Time [s]')
ylabel('Valve Opening')

fname=[fname_base,'ur'];
    
if exist('saveplots','var') && (saveplots ~= 0)
    saveas(figs{1},[fname,'.fig']);
    saveas(figs{1},[fname,'.pdf']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
    fname,'.tex'],'width','0.8\linewidth','figurehandle',figs{1},'showInfo',false);
end

figs = figure;
plot(t_dist,u_dist)
title('Tank Output Valve Opening')
ylabel('Valve Opening')
xlabel('Time [s]')
ylim([0, max(u_dist)*1.1])

%% Cost function

figs = PlotResults(results,'J');
title('Cost Function')
xlabel('Time [s]')
ylabel('Cost Function')

fname=[fname_base,'j'];
    
if exist('saveplots','var') && (saveplots ~= 0)
    saveas(figs{1},[fname,'.fig']);
    saveas(figs{1},[fname,'.pdf']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
    fname,'.tex'],'width','0.8\linewidth','figurehandle',figs{1},'showInfo',false);
end


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

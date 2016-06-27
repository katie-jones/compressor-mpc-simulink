filename = '/home/katie/school/MasterThesis/cpp/parallel/cent_output.dat';
fname_info = '/home/katie/school/MasterThesis/cpp/parallel/cent_info.dat';
res_cent = read_cpp_results(filename,fname_info,11,4,4);

n_solver_iterations = 3;

filename = ['/home/katie/school/MasterThesis/cpp/parallel/coop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/parallel/coop_info',num2str(n_solver_iterations),'.dat'];
res_coop = read_cpp_results(filename,fname_info,11,4,4);

filename = ['/home/katie/school/MasterThesis/cpp/parallel/noncoop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/parallel/noncoop_info',num2str(n_solver_iterations),'.dat'];
res_ncoop = read_cpp_results(filename,fname_info,11,4,4);

%% Cost function
res_cent = AddCostFunction(res_cent);
res_coop = AddCostFunction(res_coop);
res_ncoop = AddCostFunction(res_ncoop);

results = {res_cent, res_coop, res_ncoop};

inds = res_cent.t>=50.0;
xlims = [0 300];

norm_sd = sum(res_cent.Jsd);
norm_p = sum(res_cent.Jp);

fig=PlotCostFunction(results,'Jsd','Surge Distance');
set(fig,'units','normalized','position',[0,0.5,0.3,0.4])
legend('Centralized','Cooperative (3 it.)','Non-cooperative (3 it.)')
xlim(xlims);

fig=PlotCostFunction(results,'Jp','Output Pressure');
set(fig,'units','normalized','position',[0.3,0.5,0.3,0.4])
xlim(xlims);

fig=PlotCostFunction(results,'Jtd','Torque Setting');
set(fig,'units','normalized','position',[0,0,0.3,0.4])
xlim(xlims);

fig=PlotCostFunction(results,'Jur','Recycle Valve Opening');
set(fig,'units','normalized','position',[0.3,0,0.3,0.4])
xlim(xlims);

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



%%
figure; 
plot(res_cent.t,[res_cent.y(:,1),res_coop.y(:,1),res_ncoop.y(:,1)]);
legend('Centralized','Cooperative','Non-cooperative')
title('Surge Distance Response')

%%
figure; 
plot(res_cent.t,[res_cent.y(:,4),res_coop.y(:,4),res_ncoop.y(:,4)]);
legend('Centralized','Cooperative','Non-cooperative')
title('Output Pressure Response')

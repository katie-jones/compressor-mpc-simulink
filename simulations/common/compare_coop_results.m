N = 9;
res_coop = cell(N,1);
res_ncoop = cell(N,1);

legstr = cell(N,1);

for n_solver_iterations=1:N
    
filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/coop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/parallel/output/coop_info',num2str(n_solver_iterations),'.dat'];
res = read_cpp_results(filename,fname_info,11,4,4);
res_coop{n_solver_iterations} = AddCostFunction(res);

filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/noncoop_output',num2str(n_solver_iterations),'.dat'];
fname_info = ['/home/katie/school/MasterThesis/cpp/parallel/output/noncoop_info',num2str(n_solver_iterations),'.dat'];
res = read_cpp_results(filename,fname_info,11,4,4);
res_ncoop{n_solver_iterations} = AddCostFunction(res);

legstr{n_solver_iterations} = [num2str(n_solver_iterations),' solver it.'];
end
res_coop_orig = res_coop;
legstr_orig = legstr;

%% Cost function
inds = res_coop{1}.t>=50.0;
xlims = [50 300];

res_coop = res_coop_orig(1:2:end);
legstr = legstr_orig(1:2:end);

fig=PlotCostFunction(res_coop,'Jsd','Surge Distance',0);
set(fig,'units','normalized','position',[0,0.5,0.3,0.4])
xlim(xlims);
legend(legstr);
if saveplots
    saveas(fig,'sd_coop_penalty.fig');
    saveas(fig,'sd_coop_penalty.pdf');
end

fig=PlotCostFunction(res_coop,'Jp','Output Pressure',0);
set(fig,'units','normalized','position',[0.3,0.5,0.3,0.4])
xlim(xlims);
if saveplots
    saveas(fig,'pout_coop_penalty.fig');
    saveas(fig,'pout_coop_penalty.pdf');
end

fig=PlotCostFunction(res_coop,'Jtd','Torque Setting',0);
set(fig,'units','normalized','position',[0,0,0.3,0.4])
xlim(xlims);
if saveplots
    saveas(fig,'td_coop_penalty.fig');
    saveas(fig,'td_coop_penalty.pdf');
end

fig=PlotCostFunction(res_coop,'Jur','Recycle Valve Opening',0);
set(fig,'units','normalized','position',[0.3,0,0.3,0.4])
xlim(xlims);
if saveplots
    saveas(fig,'ur_coop_penalty.fig');
    saveas(fig,'ur_coop_penalty.pdf');
end


%%
% figure;
% plot(res_coop{1}.t,[res_coop{1}.y(:,1),res_coop{3}.y(:,1),res_coop{9}.y(:,1)]);
% legend('Centralized','Cooperative','Non-cooperative')
% title('Surge Distance Response')

%%
% figure; 
% set(gca,'YScale','log')
% hold on
% meanval=zeros(N-1,4);
% for i=1:N-1
%     semilogy(res_coop{i}.t(inds)-50,-res_coop{i}.Jur(inds)+res_coop{N}.Jur(inds));
%     meanval(i,4) = sum(res_coop{i}.Jur) - sum(res_coop{N}.Jur);
%     meanval(i,3) = sum(res_coop{i}.Jtd) - sum(res_coop{N}.Jtd);
%     meanval(i,2) = sum(res_coop{i}.Jsd) - sum(res_coop{N}.Jsd);
%     meanval(i,1) = sum(res_coop{i}.Jp) - sum(res_coop{N}.Jp);
% end
% xlim([0 300])
% xlabel('Time [s]')
% ylabel('Difference in penalty term')
% title('Recycle Valve Opening')
% legend('1','3','5','7')
% 

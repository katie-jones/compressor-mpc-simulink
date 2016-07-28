Ts = 0.05; % sampling time
Tsim = 500; % simulation time
N = 10; % max number of iterations
n_runs = 1; % number of trials

% coop_times = zeros(500/Ts,N);
coop_times = zeros(1,N);
noncoop_times = coop_times;

% cent_times = zeros(500/Ts,1);
cent_times = 0;
cent_0_times = cent_times;

for j=1:n_runs
    filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/cent_cpu_times_run',num2str(j),'.dat'];
    f = fopen(filename,'r');
    A = fscanf(f,'%g'); 
    fclose(f);
    cent_times = cent_times + A;%A(2:end)-A(1:end-1);

    filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/cent_cpu_times0_run',num2str(j),'.dat'];
    f = fopen(filename,'r');
    A = fscanf(f,'%g'); 
    fclose(f);
    cent_0_times = cent_0_times + A;%A(2:end)-A(1:end-1);
    
    for i=1:N
        filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/coop_cpu_times',num2str(i-1),'_run',num2str(j),'.dat'];
        f = fopen(filename,'r');
        A = fscanf(f,'%g');
        fclose(f);
        coop_times(:,i) = coop_times(:,i) + A;%A(2:end)-A(1:end-1);
        
        filename = ['/home/katie/school/MasterThesis/cpp/parallel/output/ncoop_cpu_times',num2str(i-1),'_run',num2str(j),'.dat'];
        f = fopen(filename,'r');
        A = fscanf(f,'%g');
        fclose(f);
        noncoop_times(:,i) = noncoop_times(:,i) + A;%A(2:end)-A(1:end-1);
    end
end

coop_times = coop_times./n_runs;
noncoop_times = noncoop_times./n_runs;
cent_times = cent_times./n_runs;
cent_0_times = cent_0_times./n_runs;

cent_mean_times = [mean(cent_0_times); mean(cent_times)*ones(N-1,1)]/Tsim * Ts;
noncoop_mean_times = noncoop_times/2/Tsim*Ts;
coop_mean_times = coop_times/2/Tsim*Ts;
% noncoop_mean_times = mean(noncoop_times);
% coop_mean_times = mean(coop_times);

%%
fig=figure; 
plot(1:N-1,cent_mean_times(2:end)/1.0e6,'-o', 1:N-1, coop_mean_times(2:end)/1e6,...
    '-s',1:N-1,noncoop_mean_times(2:end)/1e6,'-x','linewidth',1.5,'markersize',10);
grid on
% title('Controller Computation Time per Iteration')
title('Parallel System')
xlabel('Number of solver iterations')
ylabel('Average computation time [ms]')
legend('Centralized','Cooperative','Non-cooperative','location','northwest')
fig=printplot(fig);
fname='parallel_computation_cost';
if exist('saveplots','var') && saveplots~=0
    saveas(fig,[fname,'.fig'])
    saveas(fig,[fname,'.pdf'])
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
        fname,'.tex'],'width','0.8\linewidth','figurehandle',fig,'showInfo',false);
end

fig=figure;
plot(1:N-1, 1 - cent_mean_times(1)./cent_mean_times(2:end), '-o',1:N-1,...
    1- coop_mean_times(1)./coop_mean_times(2:end),'-s',1:N-1,...
    1-noncoop_mean_times(1)./noncoop_mean_times(2:end),'-x','linewidth',1.5,'markersize',10);
grid on
% title({'Relative Computation Time of QP Solver','As % of Total'})
title('Parallel System')
xlabel('Number of solver iterations')
ylabel('Average computation time [% of total]')
legend('Centralized','Cooperative','Non-cooperative','location','northwest')
fig=printplot(fig);
fname = 'parallel_qp_cost';
if exist('saveplots','var') && saveplots~=0
    saveas(fig,[fname,'.fig'])
    saveas(fig,[fname,'.pdf'])
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
        fname,'.tex'],'width','0.8\linewidth','figurehandle',fig,'showInfo',false);
end


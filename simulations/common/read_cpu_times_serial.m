Ts = 0.05; % sampling time
N = 9; % max number of iterations

coop_times = zeros(500/Ts,N);
coop_smoothed_times = zeros(floor(size(coop_times,1)/10),N);

noncoop_times = coop_times;
noncoop_smoothed_times = coop_smoothed_times;

for i=1:N
    filename = ['/home/katie/school/MasterThesis/cpp/serial/output/coop_cpu_times',num2str(i),'.dat'];
    f = fopen(filename,'r');
    A = fscanf(f,'%g'); 
    fclose(f);
    coop_times(:,i) = A(2:end)-A(1:end-1);
    for j=1:size(coop_smoothed_times,1)
        coop_smoothed_times(j,i) = mean(coop_times((j-1)*10+1:j*10,i));
    end
    
    filename = ['/home/katie/school/MasterThesis/cpp/serial/output/ncoop_cpu_times',num2str(i),'.dat'];
    f = fopen(filename,'r');
    A = fscanf(f,'%g'); 
    fclose(f);
    noncoop_times(:,i) = A(2:end)-A(1:end-1);
    for j=1:size(coop_smoothed_times,1)
        noncoop_smoothed_times(j,i) = mean(coop_times((j-1)*10+1:j*10,i));
    end
end

filename = '/home/katie/school/MasterThesis/cpp/serial/output/cent_cpu_times.dat';
f = fopen(filename,'r');
A = fscanf(f,'%g'); 
fclose(f);
cent_times = A(2:end)-A(1:end-1);

cent_mean_times = mean(cent_times)*ones(N,1);
noncoop_mean_times = mean(noncoop_times);
coop_mean_times = mean(coop_times);

%%
fig=figure; 
plot(1:N,[cent_mean_times, coop_mean_times(:), noncoop_mean_times(:)]/1.0e6,'-o','linewidth',1.5);
grid on
title('Comparison of Computational Cost')
xlabel('Number of solver iterations')
ylabel('Average computation time [ms]')
legend('Centralized','Cooperative','Non-cooperative','location','best')
fig=printplot(fig);
if exist('saveplots','var') && saveplots~=0
    saveas(fig,'serial_computation_cost.fig')
    saveas(fig,'serial_computation_cost.pdf')
end
cen_filename = 'results/timing/cent_cpu_times.dat';
coop_filename_root = 'results/timing/coop_cpu_times';
ncoop_filename_root = 'results/timing/noncoop_cpu_times';

f = fopen(cen_filename,'r');
t_cen_full = fscanf(f,'%f');
fclose(f);

N = 9;
Ts = 0.05;

t_coop_full = zeros(length(t_cen_full),N);
t_ncoop_full = t_coop_full;

for i=1:N
    f = fopen([coop_filename_root, num2str(i), '.dat'],'r');
    t_coop_full(:,i) = fscanf(f, '%f');
    fclose(f);
    f = fopen([ncoop_filename_root, num2str(i), '.dat'],'r');
    t_ncoop_full(:,i) = fscanf(f, '%f');
    fclose(f);
end

%%
to_skip = 30;
t_cen = t_cen_full(1:to_skip:end);
t_coop = t_coop_full(1:to_skip:end,:);
t_ncoop = t_ncoop_full(1:to_skip:end,:);

t_coop_mean = mean(t_coop,1);
t_ncoop_mean = mean(t_ncoop,1);
t_cen_mean = mean(t_cen)*ones(size(t_coop_mean));

t_sim = 0:Ts:(length(t_cen)-1)*Ts;
t_sim = t_sim * to_skip;
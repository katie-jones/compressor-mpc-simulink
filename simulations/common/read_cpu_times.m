Ts = 0.05; % sampling time
N = 9; % max number of iterations

coop_times = zeros(1000/Ts+1,N);
coop_smoothed_times = zeros(floor(size(coop_times,1)/10),N);

noncoop_times = coop_times;
noncoop_smoothed_times = coop_smoothed_times;

for i=1:N
    filename = ['/home/katie/school/MasterThesis/cpp/build/coop_cpu_times',num2str(i),'.dat'];
    f = fopen(filename,'r');
    A = fscanf(f,'%g'); 
    fclose(f);
    coop_times(:,i) = A(2:end)-A(1:end-1);
    for j=1:size(coop_smoothed_times,1)
        coop_smoothed_times(j,i) = mean(coop_times((j-1)*10+1:j*10,i));
    end
    
    filename = ['/home/katie/school/MasterThesis/cpp/build/noncoop_cpu_times',num2str(i),'.dat'];
    f = fopen(filename,'r');
    A = fscanf(f,'%g'); 
    fclose(f);
    noncoop_times(:,i) = A(2:end)-A(1:end-1);
    for j=1:size(coop_smoothed_times,1)
        noncoop_smoothed_times(j,i) = mean(coop_times((j-1)*10+1:j*10,i));
    end
end

filename = '/home/katie/school/MasterThesis/cpp/build/cent_cpu_times.dat';
f = fopen(filename,'r');
A = fscanf(f,'%g'); 
fclose(f);
cent_times = A(2:end)-A(1:end-1);
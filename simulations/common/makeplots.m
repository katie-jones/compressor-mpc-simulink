% set default plotting params
set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',1.5)

% if variable hasn't been set
if ~exist('saveplots','var')
    saveplots = 0;
end

% check disturbances
udist = [];
ulegstr = {};
for i=1:length(udist1)
    if udist1(i)~=0
        switch i
            case 1
                udist = cat(1,udist,uoff1(2)+udist1(1)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Input valve 1';
            case 2
                udist = cat(1,udist,uoff1(3)+udist1(2)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve 1';
            case 3
                udist = cat(1,udist,ud+udist1(3)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve (tank)';
            case 4
                udist = cat(1,udist,uoff2(2)+udist1(4)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Input valve 2';
            case 5
                udist = cat(1,udist,uoff2(3)+udist1(5)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve 2';
        end
        if udist2(i)~=0
            udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
        end
    elseif udist2(i)~=0
        switch i
            case 1
                udist = cat(1,udist,uoff1(2)+udist2(1)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Input valve 1';
            case 2
                udist = cat(1,udist,uoff1(3)+udist2(2)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve 1';
            case 3
                udist = cat(1,udist,ud+udist2(3)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve (tank)';
            case 4
                udist = cat(1,udist,uoff2(2)+udist2(4)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Input valve 2';
            case 5
                udist = cat(1,udist,uoff2(3)+udist2(5)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve 2';
        end
    end
end
    

fig = figure;
set(fig,'units','normalized','position',[0 0 1 1])

subplot(5,2,1)
plot(Td.time,Td.signals.values(:,1))
grid on
title('Compressor 1')
ylabel('Torque setting')

subplot(5,2,2)
plot(Td.time,Td.signals.values(:,2))
grid on
title('Compressor 2')

subplot(5,2,3)
plot(u_rec.time,u_rec.signals.values(:,1))
grid on
ylabel('Recycle opening')

subplot(5,2,4)
plot(u_rec.time,u_rec.signals.values(:,2))
grid on

subplot(5,2,5)
plot(pout.time,pout.signals.values(:,1))
grid on
ylabel('Output pressure')

subplot(5,2,6)
plot(pout.time,pout.signals.values(:,2))
grid on

subplot(5,2,7)
plot(SD.time,SD.signals.values(:,1))
grid on
hold on
plot([0 max(SD.time)], [yref(1) yref(1)],'-.k')
ylabel('Surge Distance')

subplot(5,2,8)
plot(SD.time,SD.signals.values(:,2))
grid on
hold on
plot([0 max(SD.time)], [yref(2) yref(2)],'-.k')

subplot(5,2,9)
plot([0 tdist(1) tdist(1) tdist(2) tdist(2) SD.time(end)],udist);
legend(ulegstr,'fontsize',12,'location','north','orientation','horizontal');
grid on
title('Valve disturbance')
ylabel('Valve setting')
xlabel('Time [s]')

subplot(5,2,10)
plot(PD.time,PD.signals.values)
grid on
hold on
plot([0 max(PD.time)], [yref(4) yref(4)],'-.k')
title('Output tank pressure')
ylabel('Pressure')
xlabel('Time [s]')

if saveplots
    set(fig,'units','centimeters','position',1.25*[0 0 21 29.7])
    fig=printplot(fig);
    
    if (results_folder(end)~='/')
        results_folder = strcat(results_folder,'/');
    end
    if dist_dirname(end)~='/'
        dist_dirname = strcat(dist_dirname,'/');
    end
    basename = [results_folder,dist_dirname,results_fname];
    if ~exist(results_folder,'dir')
        mkdir(results_folder)
    elseif ~exist([results_folder,dist_dirname],'dir')
        mkdir([results_folder,dist_dirname])
    elseif results_overwrite==0
        
        n = 0;
        while (exist([basename,'.pdf'],'file') || exist([basename,'.mat'],'file') || exist([basename,'.fig'],'file'))
            n=n+1;
            basename=[basename,num2str(n)];
        end
    end
    saveas(fig,[basename,'.fig'])
    saveas(fig,[basename,'.pdf'])
    [n_delay,~,~,p,m,UWT,YWT] = const_mpc();
    [n_barrier,delta_barrier] = const_barrier();
    Results.n_delay = n_delay;
    Results.p = p;
    Results.m = m;
    Results.UWT = UWT;
    Results.YWT = YWT;
    Results.lb = lb;
    Results.ub = ub;
    Results.Td = Td;
    Results.ur = u_rec;
    Results.SD = SD;
    Results.pout = pout;
    Results.pd = PD;
    Results.tdist = tdist;
    Results.udist = [udist1; udist2];
    Results.n_barrier = n_barrier;
    Results.delta_barrier = delta_barrier;
    save([basename,'.mat'],'Results')
end
    
set(0,'defaultlinelinewidth',1)
    

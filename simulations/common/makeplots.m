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
input_added = [0 0];
output_added = [0 0];
for i=1:length(udist1)
    if udist1(i)~=0
        switch mod(i,3)
        case 1
            if ~input_added(1)
                input_added(1) = 1;
                udist = cat(1,udist,uoff1(2)+udist1(1)*[0 0 1 1 1 1]);
                if udist2(1)~=0
                    udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
                end
                ulegstr{end+1} = 'Input valve 1';
                udist = cat(1,udist,uoff2(2)+udist1(4)*[0 0 1 1 1 1]);
                if udist2(4)~=0
                    udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
                end
                ulegstr{end+1} = 'Input valve 2';
            end
        case 2
            if ~output_added(1)
                output_added(1) = 1;
                udist = cat(1,udist,uoff1(3)+udist1(2)*[0 0 1 1 1 1]);
                if udist2(2)~=0
                    udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
                end
                ulegstr{end+1} = 'Output valve 1';
                udist = cat(1,udist,uoff2(3)+udist1(5)*[0 0 1 1 1 1]);
                if udist2(5)~=0
                    udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
                end
                ulegstr{end+1} = 'Output valve 2';
            end
        case 0
            udist = cat(1,udist,ud+udist1(3)*[0 0 1 1 1 1]);
            if udist2(3)~=0
                udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
            end
            ulegstr{end+1} = 'Output valve (tank)';
        end
    elseif udist2(i)~=0
        switch mod(i,3)
        case 1
            if ~input_added(2)
                input_added(2) = 1;
                udist = cat(1,udist,uoff1(2)+udist2(1)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Input valve 1';
                udist = cat(1,udist,uoff2(2)+udist2(4)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Input valve 2';
            end
        case 2
            if ~output_added(2)
                output_added(2) = 1;
                udist = cat(1,udist,uoff1(3)+udist2(2)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve 1';
                udist = cat(1,udist,uoff2(3)+udist2(5)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve 2';
            end
        case 3
            udist = cat(1,udist,ud+udist2(3)*[0 0 1 1 1 1]);
            ulegstr{end+1} = 'Output valve (tank)';
        end
    end
end


fig = figure;
if saveplots
    set(fig,'units','centimeters','position',1.25*[0 0 21 29.7])
    fig=printplot(fig);
else
    set(fig,'units','normalized','position',[0 0 1 1])
end

subplot(3,2,1)
plot(Td.time,Td.signals.values)
grid on
title('Inputs','fontsize',16)
ylabel('Torque setting')
legend('Comp. 1','Comp. 2')

subplot(3,2,2)
plot(pout.time,pout.signals.values)
grid on
title('Outputs','fontsize',16)
ylabel('Output pressure')
legend('Comp. 1','Comp. 2')

subplot(3,2,3)
plot(u_rec.time,u_rec.signals.values)
grid on
ylabel('Recycle opening')
legend('Comp. 1','Comp. 2','location','best')

subplot(3,2,4)
plot(SD.time,SD.signals.values)
grid on
hold on
plot([0 max(SD.time)], [yref(1) yref(1)],'-.k')
ylabel('Surge Distance')
legend('Comp. 1','Comp. 2','location','southeast')

subplot(3,2,5)
plot([0 tdist(1) tdist(1) tdist(2) tdist(2) SD.time(end)],udist);
legend(ulegstr,'fontsize',12,'location','north','orientation','horizontal');
grid on
ylabel('Valve setting')
xlabel('Time [s]')

subplot(3,2,6)
plot(PD.time,PD.signals.values)
grid on
hold on
plot([0 max(PD.time)], [yref(4) yref(4)],'-.k')
ylabel('Pressure')
xlabel('Time [s]')


[n_delay,~,~,p,m] = const_mpc();
weights;
[n_barrier,delta_barrier] = const_barrier();
[~,Vtank] = const_tank();
Results.n_delay = n_delay;
Results.Vtank = Vtank;
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
Results.yref = yref;
Results.n_barrier = n_barrier;
Results.delta_barrier = delta_barrier;
if exist('n_iterations','var')
Results.n_iterations = n_iterations;
end
Results.label = results_label;

if saveplots


    if (results_folder(end)~='/')
        results_folder = strcat(results_folder,'/');
    end
    if dist_dirname(end)~='/'
        dist_dirname = strcat(dist_dirname,'/');
    end
    basename = [results_folder,dist_dirname,results_fname];
    if ~exist(results_folder,'dir')
        mkdir(results_folder)
    end
    if ~exist([results_folder,dist_dirname],'dir')
        mkdir([results_folder,dist_dirname])
    elseif results_overwrite==0


        if (exist([basename,'.pdf'],'file') || exist([basename,'.mat'],'file') || exist([basename,'.fig'],'file'))
            n = 1;
            while (exist([basename,num2str(n),'.pdf'],'file') || exist([basename,num2str(n),'.mat'],'file') || exist([basename,num2str(n),'.fig'],'file'))
                n=n+1;
            end
            basename=[basename,num2str(n)];

        end
    end
    saveas(fig,[basename,'.fig'])
    saveas(fig,[basename,'.pdf'])

    save([basename,'.mat'],'Results')
    print_results_info(results_folder,results_label);
end

set(0,'defaultlinelinewidth',1)


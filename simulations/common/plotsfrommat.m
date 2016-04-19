% set default plotting params
set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',1.5)

% if variable hasn't been set
if ~exist('saveplots','var')
    saveplots = 0;
end

v2struct(Results)

udist1 = udist(1,:);
udist2 = udist(2,:);

if ~exist('off1','var')
    [~, ~, ~, ~, ~, uoff1, uoff2, ud] = const_sim();
end

if ~exist('yref','var')
    yref = [0.2175 0.2175 0 1.12]';
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
                    ulegstr{end+1} = 'Input valve 1';
                    udist = cat(1,udist,uoff2(2)+udist1(4)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Input valve 2';
                end
            case 2
                if ~output_added(1)
                    output_added(1) = 1;
                    udist = cat(1,udist,uoff1(3)+udist1(2)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Output valve 1';
                    udist = cat(1,udist,uoff2(3)+udist1(5)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Output valve 2';
                end
            case 3
                udist = cat(1,udist,ud+udist1(3)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve (tank)';
        end
        if udist2(i)~=0
            udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
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
set(fig,'units','normalized','position',[0 0 1 1])

subplot(3,2,1)
plot(Td.time,Td.signals.values)
grid on
title('Inputs','fontsize',16)
ylabel('Torque Setting')
legend('Comp. 1','Comp. 2','location','best')

subplot(3,2,2)
plot(pout.time, pout.signals.values)
grid on
title('Outputs','fontsize',16)
ylabel('Output pressure')
legend('Comp. 1','Comp. 2','location','best')

subplot(3,2,3)
plot(ur.time, ur.signals.values)
grid on
ylabel('Recycle valve')
legend('Comp. 1','Comp. 2','location','best')

subplot(3,2,4)
plot(SD.time, SD.signals.values)
grid on
hold on
plot([0 max(SD.time)], [yref(2) yref(2)],'-.k')
ylabel('Surge distance')
legend('Comp. 1','Comp. 2','Ref','location','best')

subplot(3,2,5)
plot([0 tdist(1) tdist(1) tdist(2) tdist(2) SD.time(end)],udist);
grid on
ylabel('Valve setting')
legend(ulegstr,'location','best')

subplot(3,2,6)
plot(pd.time,pd.signals.values)
grid on
ylabel('Tank pressure')
hold on
plot([0 max(pd.time)], [yref(4) yref(4)],'-.k')


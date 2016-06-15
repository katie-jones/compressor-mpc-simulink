function [] = plotcompare(res1, res2, res3, tlims)
% set default plotting params
set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',1.5)

if sum(sum(res1.udist~=res2.udist))~=0
    error('Two datasets have different disturbances')
elseif sum(sum(res1.tdist~=res2.tdist))~=0
    error('Two datasets have different disturbance times')
end

if nargin < 4
    tlims = [0, max(res1.SD.time)];
end

tdist = res1.tdist;

% if symmetric, only plot one response
    inds = 1:2;
    legstr = {res1.label, res2.label};

udist1 = res1.udist(1,:);
udist2 = res1.udist(2,:);

if ~exist('res1.off1','var')
    [~, ~, ~, ~, ~, uoff1, uoff2] = const_sim();
end

yref = res1.yref;

% check disturbances
udist = [];
ulegstr = {};
input_added = [0 0];
output_added = [0 0];
for i=1:length(udist1)
    if udist1(i)~=0
        switch i
            case 0
                    input_added(1) = 1;
                    udist = cat(1,udist,uoff1(2)+udist1(1)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Input valve 1';
                    if udist2(1)~=0
                        udist(end,:) = udist(end,:) + udist2(i)*[0 0 0 0 1 1];
                    end
                
            case 1
                
                    output_added(1) = 1;
                    udist = cat(1,udist,uoff1(3)+udist1(2)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Output valve 1';
                    if udist2(2)~=0
                        udist(end,:) = udist(end,:) + udist2(2)*[0 0 0 0 1 1];
                    end
            case 2
                    udist = cat(1,udist,uoff2(3)+udist1(3)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Output valve 2';
                    if udist2(3)~=0
                        udist(end,:) = udist(end,:) + udist2(3)*[0 0 0 0 1 1];
                    end
                
        end
    elseif udist2(i)~=0
        switch i
            case 0
                
                    
                    udist = cat(1,udist,uoff1(2)+udist2(1)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Input valve 1';
                
            case 1
                
                    
                    udist = cat(1,udist,uoff1(3)+udist2(2)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Output valve 1';
            case 2
                    udist = cat(1,udist,uoff2(3)+udist2(3)*[0 0 1 1 1 1]);
                    ulegstr{end+1} = 'Output valve 2';
                
        end
    end
end


fig = figure;
set(fig, 'units','centimeters','position',[0 0 25 20]);

subplot(2,2,1)
plot(res1.Td.time, res1.Td.signals.values(:,1), res2.Td.time, res2.Td.signals.values(:,1))
grid on
title('Compressor 1','fontsize',16)
ylabel('Torque Setting')
legend(legstr,'location','best')
xlim([tlims])

subplot(2,2,2)
plot(res1.Td.time, res1.Td.signals.values(:,2), res2.Td.time, res2.Td.signals.values(:,2))
grid on
title('Compressor 2','fontsize',16)
ylabel('Torque Setting')
legend(legstr,'location','best')
xlim([tlims])

subplot(2,2,3)
plot(res1.ur.time, res1.ur.signals.values(:,1),res2.ur.time,res2.ur.signals.values(:,1))
grid on
ylabel('Recycle valve')
legend(legstr,'location','best')
xlim([tlims])

subplot(2,2,4)
plot(res1.ur.time, res1.ur.signals.values(:,2),res2.ur.time,res2.ur.signals.values(:,2))
grid on
ylabel('Recycle valve')
legend(legstr,'location','best')
xlim([tlims])

fig = figure;
set(fig, 'units','centimeters','position',[0 0 25 20]);

subplot(2,2,1)
plot(res1.SD.time, res1.SD.signals.values(:,1),res2.SD.time,res2.SD.signals.values(:,1))
title('Compressor 1','fontsize',16)
grid on
hold on
plot([0 max(res1.SD.time)], [yref(2) yref(2)],'-.k')
ylabel('Surge distance')
legend(legstr,'location','best')
xlim([tlims])

subplot(2,2,2)
plot(res1.SD.time, res1.SD.signals.values(:,2),res2.SD.time,res2.SD.signals.values(:,2))
title('Compressor 2','fontsize',16)
grid on
hold on
plot([0 max(res1.SD.time)], [yref(4) yref(4)],'-.k')
ylabel('Surge distance')
legend(legstr,'location','best')
xlim([tlims])

subplot(2,2,3)
plot(res1.pout.time, res1.pout.signals.values(:,1),res2.pout.time,res2.pout.signals.values(:,1))
grid on
hold on
plot([0 max(res1.pout.time)], [yref(1) yref(1)],'-.k')
ylabel('Output pressure')
legend(legstr,'location','best')
xlim([tlims])

subplot(2,2,4)
plot(res1.pout.time, res1.pout.signals.values(:,2),res2.pout.time,res2.pout.signals.values(:,2))
grid on
hold on
plot([0 max(res1.pout.time)], [yref(3) yref(3)],'-.k')
ylabel('Output pressure')
legend(legstr,'location','best')
xlim([tlims])

end

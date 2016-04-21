function [] = plotcompare(res1, res2, res3, tlims)
% set default plotting params
set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',1.5)

if sum(sum(res1.udist~=res2.udist))~=0
    error('Two datasets have different disturbances')
elseif sum(sum(res1.udist~=res3.udist))~=0
    error('Two datasets have different disturbances')
elseif sum(sum(res1.tdist~=res2.tdist))~=0
    error('Two datasets have different disturbance times')
elseif sum(sum(res1.tdist~=res3.tdist))~=0
    error('Two datasets have different disturbance times')
end

if nargin < 4
    tlims = [0, max(res1.SD.time)];
end

tdist = res1.tdist;

% if symmetric, only plot one response
if sum(sum(res1.udist(:,1:2)~=res1.udist(:,4:5)))~=0
    inds = 1:2;
    legstr = {[res1.label,' comp. 1'],[res1.label,' comp. 2'],[res2.label,' comp. 1'],[res2.label,' comp. 2'],[res3.label,' comp. 1'],[res3.label,' comp. 2']};
else
    inds = 1;
    legstr = {res1.label,res2.label,res3.label};
end

udist1 = res1.udist(1,:);
udist2 = res1.udist(2,:);

if ~exist('res1.off1','var')
    [~, ~, ~, ~, ~, uoff1, uoff2, ud] = const_sim();
end

yref = res1.yref;

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
            case 0
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
            case 0
                udist = cat(1,udist,ud+udist2(3)*[0 0 1 1 1 1]);
                ulegstr{end+1} = 'Output valve (tank)';
        end
    end
end


fig = figure;
set(fig,'units','normalized','position',[0 0 1 1])

subplot(3,2,1)
plot(res1.Td.time,res1.Td.signals.values(:,inds),res2.Td.time,res2.Td.signals.values(:,inds),res3.Td.time,res3.Td.signals.values(:,inds))
grid on
title('Inputs','fontsize',16)
ylabel('Torque Setting')
legend(legstr,'location','best')
xlim([tlims])

subplot(3,2,2)
plot(res1.pout.time, res1.pout.signals.values(:,inds),res2.pout.time,res2.pout.signals.values(:,inds),res3.pout.time,res3.pout.signals.values(:,inds))
grid on
title('Outputs','fontsize',16)
ylabel('Output pressure')
legend(legstr,'location','best')
xlim([tlims])

subplot(3,2,3)
plot(res1.ur.time, res1.ur.signals.values(:,inds),res2.ur.time,res2.ur.signals.values(:,inds),res3.ur.time,res3.ur.signals.values(:,inds))
grid on
ylabel('Recycle valve')
legend(legstr,'location','best')
xlim([tlims])

subplot(3,2,4)
plot(res1.SD.time, res1.SD.signals.values(:,inds),res2.SD.time,res2.SD.signals.values(:,inds),res3.SD.time,res3.SD.signals.values(:,inds))
grid on
hold on
plot([0 max(res1.SD.time)], [yref(2) yref(2)],'-.k')
ylabel('Surge distance')
legend(legstr,'location','best')
xlim([tlims])

subplot(3,2,5)
plot([0 tdist(1) tdist(1) tdist(2) tdist(2) res1.SD.time(end)],udist);
grid on
ylabel('Valve setting')
legend(ulegstr,'location','best')
xlim([tlims])

subplot(3,2,6)
plot(res1.pd.time,res1.pd.signals.values,res2.pd.time,res2.pd.signals.values,res3.pd.time,res3.pd.signals.values)
grid on
ylabel('Tank pressure')
hold on
plot([0 max(res1.pd.time)], [yref(4) yref(4)],'-.k')
legend(res1.label,res2.label,res3.label)
xlim([tlims])

end

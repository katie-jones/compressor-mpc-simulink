%% Get data
n_max_iter = 20;

folder='/home/katie/school/MasterThesis/cpp/serial/output/';
f = fopen([folder,'coop_J',num2str(n_max_iter),'.dat'],'r');
A = fscanf(f,'%f',[4*n_max_iter,Inf]);
fclose(f);

Ju = A(1:2:end,:)';
Jy = A(2:2:end,:)';

Ju1 = Ju(:,1:2:end);
Ju2 = Ju(:,2:2:end);

Jy1 = Jy(:,1:2:end);
Jy2 = Jy(:,2:2:end);

J1 = 0.5*Ju1 + Jy1;
J2 = 0.5*Ju2 + Jy2;

%% Make plots

fname_base = 'serial_';

j1_zoom = J1(1000:1300,:);
j2_zoom = J2(1000:1300,:);

n = 0:length(j1_zoom)-1;
Ts = 0.05;
t = 50 + n*Ts;

j1min = min(j1_zoom,[],2);

dj1_zoom = (j1_zoom - repmat(j1_zoom(:,1),1,n_max_iter))./max(abs(repmat(j1_zoom(:,end),1,n_max_iter)),0.5)*100;
dj2_zoom = (j2_zoom - repmat(j2_zoom(:,1),1,n_max_iter))./max(abs(repmat(j2_zoom(:,end),1,n_max_iter)),0.5)*100;

figure;
hold on; grid on
plot(mean(dj1_zoom))
plot(mean(dj2_zoom))
set(gca,'ColorOrderIndex',1)
plot(3,mean(dj1_zoom(:,3)),'+k','markersize',14)
plot(3,mean(dj2_zoom(:,3)),'+k','markersize',14)
xlabel('No. Solver Iterations')
ylabel('Reduction in cost function [%]')
title({'Average Normalized Reduction in Cost Function','Serial Cooperative Controller'})
xlim([0 15])

legend('Compressor 1','Compressor 2','Operating Point')

fname=[fname_base,'dj'];

if exist('saveplots','var') && (saveplots ~= 0)
    saveas(gcf,[fname,'.fig']);
    saveas(gcf,[fname,'.pdf']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
    fname,'.tex'],'width','0.45\linewidth','showInfo',false);
end

% 
% figure; 
% semilogy(t,(dj1_zoom(:,[1,2,3,4])));
% hold on;
% plot(t(end),dj1_zoom(end,1),'-k')
% plot(t(end),dj1_zoom(end,1),'--k')
% set(gca,'ColorOrderIndex',1)
% semilogy(t,-(dj1_zoom(:,[1,2,3,4])),'--')
% legend('1 it.','2 it.','3 it.','4 it.','Positive change','Negative change','location','south','orientation','horizontal')
% title('Normalized Change in Cost Function');
% ylabel('Absolute change in cost function')
% 
% fname=[fname_base,'dj1'];
% 
% if exist('saveplots','var') && (saveplots ~= 0)
%     saveas(gcf,[fname,'.fig']);
%     saveas(gcf,[fname,'.pdf']);
%     matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
%     fname,'.tex'],'width','0.8\linewidth','showInfo',false);
% end
% 
% figure;
% plot(t,j1_zoom(:,20))
% title('Cost Function at 20 Iterations')
% ylabel('Cost Function')
% xlabel('Time [s]')
% fname=[fname_base,'j1'];
% 
% if exist('saveplots','var') && (saveplots ~= 0)
%     saveas(gcf,[fname,'.fig']);
%     saveas(gcf,[fname,'.pdf']);
%     matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
%     fname,'.tex'],'width','0.8\linewidth','showInfo',false);
% end
% 
% figure; 
% semilogy(t,(dj2_zoom(:,[1,2,3,4])));
% hold on;
% plot(t(end),dj2_zoom(end,1),'-k')
% plot(t(end),dj2_zoom(end,1),'--k')
% set(gca,'ColorOrderIndex',1)
% semilogy(t,-(dj2_zoom(:,[1,2,3,4])),'--')
% legend('1 it.','2 it.','3 it.','4 it.','Positive change','Negative change','location','south','orientation','horizontal')
% title('Normalized Change in Cost Function');
% ylabel('Absolute change in cost function')
% fname=[fname_base,'dj2'];
% 
% if exist('saveplots','var') && (saveplots ~= 0)
%     saveas(gcf,[fname,'.fig']);
%     saveas(gcf,[fname,'.pdf']);
%     matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
%     fname,'.tex'],'width','0.8\linewidth','showInfo',false);
% end
% 
% figure;
% plot(t,j2_zoom(:,20))
% title('Cost Function at 20 Iterations')
% ylabel('Cost Function')
% xlabel('Time [s]')
% fname=[fname_base,'j2'];
% 
% if exist('saveplots','var') && (saveplots ~= 0)
%     saveas(gcf,[fname,'.fig']);
%     saveas(gcf,[fname,'.pdf']);
%     matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
%     fname,'.tex'],'width','0.8\linewidth','showInfo',false);
% end

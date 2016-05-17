% [fig,map] = imread('compmap.png','png');
% imshow(fig,map);

flow_cal = [0 0.6]';
pratio_cal = [1.05 1.4]';

flow_plot_vals = linspace(-0.05,flow_cal(2),100);

if ~exist('flow_cal_coords','var')
    disp('Click on points x=0 and x=0.6');
    [flow_cal_coords,~] = ginput(2);
end

if ~exist('pratio_cal_coords','var')
    disp('Click on points y=1.05 and y=1.4');
    [~,pratio_cal_coords] = ginput(2);
end

if ~exist('surge_flow','var')
    disp('Click on surge line')
    [surge_flow,surge_pratio] = ginput_custom(100,'crosshair');
end

speeds = [1675.5161, 1884.9556, 2094.3951, 2303.8346, 2617.9939, 2932.1531];

coeffs_flow = polyfit(flow_cal_coords, flow_cal,1);
coeffs_pratio = polyfit(pratio_cal_coords, pratio_cal,1);

pratio_plot_vals = zeros(length(flow_plot_vals),length(speeds));

if ~exist('flows','var')
    
    for i=1:length(speeds)
        fprintf('Click on points from speed=%5.3g\n',speeds(i));
        [flow_in, pratio_in] = ginput_custom(100,'crosshair');
        flows{i} = flow_in;
        pratios{i} = pratio_in;
    end
    save('compmap.mat','flow_cal','pratio_cal','flow_cal_coords','pratio_cal_coords','surge_flow','surge_pratio','flows','pratios')
end

for i=1:length(speeds)
    pratio_plot_vals(:,i) = interp1(polyval(coeffs_flow,flows{i}),polyval(coeffs_pratio,pratios{i}),flow_plot_vals);
end


%%
set(0,'defaultlinelinewidth',1.5)
set(0, 'defaultlinemarkersize',6)
fig=figure;
set(fig,'units','centimeter','position',[5 5 20 20])
ax=axes; hold on;
set(ax,'fontsize',13)
plot(ax,flow_plot_vals,pratio_plot_vals,'-b')

ylims = get(ax,'ylim');
plot(ax,[0 0],[ylims(1),ylims(2)],'-k');

surge_flow_plot = [0.206,0.273];
surge_pratio_plot = interp1(polyval(coeffs_flow,surge_flow),polyval(coeffs_pratio,surge_pratio),surge_flow_plot,'linear','extrap');

plot(ax,surge_flow_plot,surge_pratio_plot,'--k', 'linewidth',2.5)

surgecontrol_flow_plot = [0.303 0.375];
plot(ax,surgecontrol_flow_plot,interp1(surge_flow_plot+0.1,surge_pratio_plot,surgecontrol_flow_plot,'linear','extrap'),'--k', 'linewidth',2.5)


for i=1:length(speeds)
    text(0.02,interp1(polyval(coeffs_flow,flows{i}),polyval(coeffs_pratio,pratios{i}),0.02)-0.01,sprintf('%.0f RPM',roundsd(speeds(i),3)),'fontsize',11,'color','blue');
end

text(0.195,1.195,'Surge line','Rotation',81,'fontsize',13,'fontweight','bold')
text(0.29,1.17,'Surge control line','Rotation',81,'fontsize',13,'fontweight','bold')


xlabel('Mass Flow Rate [kg s^{-1}]')
ylabel('Pressure Ratio')


% Plot operating points
mu_flow = 0.45;
mu_pratio = 1.35;

sig_flow = 0.025;
sig_pratio = 0.011;

N = 100;

flow_op = normrnd(mu_flow, sig_flow, N, 1);
pratio_op = normrnd(mu_pratio, sig_pratio, N, 1);

plot(flow_op, pratio_op, '+c')

set(0,'defaultlinelinewidth',1)

% %%% coefficients of the compressor map
% A = [0.000299749505193654     -0.000171254191089237      3.57321648097597e-05    -9.1783572200945e-07 ...
%     -0.252701086129365         0.136885752773673         -0.02642368327081        0.00161012740365743 ...
%     54.8046725371143          -29.9550791497765          5.27827499839098         0.693826282579158];
%
%
% mc_vec=linspace(-0.05,0.5,100);
% wc_vec = linspace(300,600,6);
%
% p_ratio = zeros(length(wc_vec),length(mc_vec));
%
% for i=1:length(wc_vec)
%     wc = wc_vec(i);
%     for j=1:length(mc_vec)
%         mc = mc_vec(j);
%
%         M = [wc^2*mc^4 wc^2*mc^3 wc^2*mc^2 wc^2*mc wc^2 ...
%             wc*mc^4 wc*mc^3 wc*mc^2 wc*mc wc ...
%             mc^4 mc^3 mc^2 mc 1]';
%         p_ratio(i,j) = A * M;
%
%     end
% end
%
% figure;
% plot(mc_vec,p_ratio)
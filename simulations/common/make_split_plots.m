function [] = make_split_plots(results)
figure;
ax1=subplot(2,2,1);
plot(results.t,results.u(:,[1,3]));
title('Torque')

ax2=subplot(2,2,2);
plot(results.t,results.u(:,[2,4]));
title('Recycle valve')

ax3=subplot(2,2,3);
plot(results.t,results.y(:,[1,3]));
title('Pressure')

ax4=subplot(2,2,4);
plot(results.t,results.y(:,[2,4]));
title('Surge Distance')

linkaxes([ax1,ax2,ax3,ax4],'x')

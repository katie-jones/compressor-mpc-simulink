load mpcrefdata
close all
figure; 
plot(u_ref.Time, u_ref.Data(:,1), u_new.Time, u_new.Data(:,1))
legend('Old controller','New controller')
title('Torque')
xlabel('Time (s)')

figure; 
plot(u_ref.Time, u_ref.Data(:,2), u_new.Time, u_new.Data(:,2))
legend('Old controller','New controller')
title('Recycle opening')
xlabel('Time (s)')

figure;
plot(y_ref.Time, y_ref.Data(:,1), y_new.Time, y_new.Data(:,1))
hold on
plot(y_ref.Time, y_ref.Data(:,3),':k','linewidth',2)
legend('Old controller','New controller','Reference signal')
title('P_{out}')
xlabel('Time (s)')

figure;
plot(y_ref.Time, y_ref.Data(:,2), y_new.Time, y_new.Data(:,2))
hold on
plot(y_ref.Time, y_ref.Data(:,4),':k','linewidth',2)
legend('Old controller','New controller','Reference signal')
title('Surge Distance')
xlabel('Time (s)')
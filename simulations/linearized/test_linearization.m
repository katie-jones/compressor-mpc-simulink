N = 200;
rec_open = 0.1*ones(N,1);
Td = linspace(-0.05,0.05,N)';
Inflow_opening = 0.405;
Outflow_opening = 0.393;

Ts = 0.05; % 50 ms

x_init_lin = [0.898
      1.126
      0.15
      439.5
      0];
u_init = [0.304,Inflow_opening,Outflow_opening,0.1]';
  
x = x_init_lin;
[Ac,Bc,Cc] = get_linearized_matrices(x, u_init);
f_initc = get_comp_deriv(x,u_init)';
[A,B,C,f_init] = discretize_rk4(Ac,Bc,Cc,f_initc,Ts);

sys2 = ss(A,B,C,0);

% t = linspace(0,10,N)';
% 
% y_linearized = lsim(sys,[Td,rec_open]',t);
% 
% figure; plot(t,y_linearized(:,1)); hold on
% plot(p_out.time, p_out.signals.values,'-g')

delta_x = rand(N,5)*0.05-0.025;
delta_x(:,2) = delta_x(:,2)/10; % sensitive to perturbations in p2

delta_u = [rand(N,1)*0.002-0.001, zeros(N,2), -0*ones(N,1)];%rand(N,1)-0.5];

for i=1:N
    f(i,:) = get_comp_deriv(x+delta_x(i,:)',u_init+delta_u(i,:)');
    f2(i,:) = f_init(:) + A*delta_x(i,:)' + B*delta_u(i,[1,4])';
end
close all
for i=1:5
    figure; 
    plot(1:N,[f(:,i),f2(:,i)])
end

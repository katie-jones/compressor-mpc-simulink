N = 200;
rec_open = 1*ones(N,1);
Td = linspace(-0.05,0.05,N)';
Inflow_opening = 0.405;
Outflow_opening = 0.393;

Ts = 0.05; % 50 ms

x_init_lin = [0.898
      1.126
      0.15
      439.5
      0];
u_init = [0.304,Inflow_opening,Outflow_opening,1]';
  
x = x_init_lin;
[Ac,Bc,Cc] = get_linearized_matrices(x, u_init([1,4])-[0.304,0]');
f_initc = get_comp_deriv(x,u_init)';% - [0;0;0;0;2*0.0263]';
[A,B,C,f_init] = discretize_rk4(Ac,Bc,Cc,f_initc,Ts);

sys2 = ss(A,B,C,0);

% t = linspace(0,10,N)';
% 
% y_linearized = lsim(sys,[Td,rec_open]',t);
% 
% figure; plot(t,y_linearized(:,1)); hold on
% plot(p_out.time, p_out.signals.values,'-g')

delta_x = (rand(N,5)-0.5).*(ones(N,1)*[0.01,0.01,0.05,10,0.05]);
% delta_x(:,2) = delta_x(:,2)/10; % sensitive to perturbations in p2
% delta_x(:,1:4) = 0;
% delta_x(:,4) = 0;

delta_u = [rand(N,1)*0.002-0.001, zeros(N,2), -1*ones(N,1)];%rand(N,1)-0.5];
% delta_u(1,:) = 0;

f = zeros(N,5);

for i=1:N
%     [Ac,Bc,Cc] = get_linearized_matrices(x+delta_x(i,:)', u_init([1,4])-[0.304;0]+delta_u(i,[1,4])');
    fc = get_comp_deriv(x+delta_x(i,:)',u_init+delta_u(i,:)');
%     [~,~,~,df] = discretize_rk4(Ac,Bc,Cc,fc,Ts);
    
%     f(i,:) = df;
    f(i,:) = fc;
    
    f2(i,:) = f_initc(:)  + Ac*delta_x(i,:)' + Bc*delta_u(i,[1,4])';
end
close all
for i=1:5
    figure; 
    plot(1:N,[f(:,i),f2(:,i)])
end

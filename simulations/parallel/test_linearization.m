% Test linearization of parallel compressor system
% Linearize system about a point, add offsets on states and inputs and 
% compare derivatives calculated using linearization to the exact values
N = 50; % number of time steps

Td = [0.3, -0.1]; % torque input
u_rec = [0.1, 0.5]; % recycle opening
u_d = 0.9; % discharge valve opening

P_D = 1.07; % initialize tank pressure
P_D_init = P_D;

Ts = 0.05;

usize = 5;
xsize = 5;

% linearization point
x_init_lin = [0.912; 1.17; 0.14; 465; 0];
[~,~,~,~,~,uoff] = const_sim();
u = [uoff; uoff; 0] + [Td(1); 0; 0; u_rec(1); 0; Td(2); 0; 0; u_rec(2); 0; u_d];
x = [x_init_lin; x_init_lin; P_D];

% initial derivatives
f1 = get_comp_deriv(x_init_lin,[u(1:usize);P_D],1);
f2 = get_comp_deriv(x_init_lin,[u(usize+1:2*usize);P_D],1);
ftank = get_tank_deriv(x,u);

% linearize system
[Ac,Bc,Cc] = linearize_tank(x, u);

% disturbances
delta_x = ones(N,1)*[repmat([0.01, 0.01, 0.05, 10, 0.05],1,2),0.02];
delta_x = (rand(2*xsize+1,N)-0.5).*delta_x';

delta_u = ones(N,1)*[repmat([0.01, 0, 0, 0.01, 0],1,2),0];
delta_u = (rand(2*usize+1,N)-0.5).*delta_u';

f = zeros(N,2*xsize+1);
fest = f;

for i=1:N
    P_D = P_D_init + delta_x(end,i);
    fe1 = get_comp_deriv(x_init_lin+delta_x(1:xsize,i),[u(1:usize)+delta_u(1:usize,i);P_D],1);
    fe2 = get_comp_deriv(x_init_lin+delta_x(xsize+1:2*xsize,i),[u(usize+1:2*usize)+delta_u(usize+1:2*usize,i);P_D],1);
    fetank = get_tank_deriv(x+delta_x(:,i),u+delta_u(:,i));
    
    f(i,:) = [fe1-f1; fe2-f2; fetank-ftank];
    
    fest(i,:) = Ac*delta_x(:,i) + Bc*delta_u([1,4,usize+1,usize+4],i);
    
end
close all

for i=1:11
    figure; 
    plot(1:N,[f(:,i),fest(:,i)])
end


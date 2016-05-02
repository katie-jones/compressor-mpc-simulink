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
[f1,m1] = get_comp_deriv(x_init_lin,[u(1:usize-1);P_D],1);
[f2,m2] = get_comp_deriv(x_init_lin,[u(usize+1:2*usize-1);P_D],1);
ftank = get_tank_deriv(P_D,[m1+m2;u_d]);

% linearize system
[Ac,Bc,Ccorig] = linearize_tank(x, u);
Cc = [Ccorig([2,4],:); Ccorig(1,:)-Ccorig(3,:); Ccorig(5,:)];

% disturbances
delta_x = ones(N,1)*[repmat([0.01, 0.01, 0.05, 10, 0.05],1,2),0.02];
delta_x = (rand(2*xsize+1,N)-0.5).*delta_x';

delta_u = ones(N,1)*[repmat([0.01, 0, 0, 0.01, 0],1,2),0];
delta_u = (rand(2*usize+1,N)-0.5).*delta_u';

f = zeros(N,2*xsize+1);
fest = f;

for i=1:N
    P_D = P_D_init + delta_x(end,i);
    [fe1,m1] = get_comp_deriv(x_init_lin+delta_x(1:xsize,i),[u(1:usize-1)+delta_u(1:usize-1,i);P_D],1);
    [fe2,m2] = get_comp_deriv(x_init_lin+delta_x(xsize+1:2*xsize,i),[u(usize+1:2*usize-1)+delta_u(usize+1:2*usize-1,i);P_D],1);
    fetank = get_tank_deriv(x(end)+delta_x(end,i),[m1+m2,u(end)+delta_u(end,i)]);
    
    f(i,:) = [fe1; fe2; fetank];
    
    fest(i,:) = [f1;f2;ftank] + Ac*delta_x(:,i) + Bc*delta_u([1,4,usize+1,usize+4],i);
    
end
close all

for i=1:11
    figure; 
    plot(1:N,[f(:,i),fest(:,i)])
end


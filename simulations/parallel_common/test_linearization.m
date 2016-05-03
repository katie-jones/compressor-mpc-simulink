% Test linearization of parallel compressor system
% Linearize system about a point, add offsets on states and inputs and 
% compare derivatives calculated using linearization to the exact values
N = 50; % number of time steps

Td = [0, 0]; % torque input
u_rec = [0, 0]; % recycle opening
u_d = 0.7; % discharge valve opening

P_D = 0.954; % initialize tank pressure
P_D_init = P_D;

[~,Pin,Pout] = const_flow();

% linearization point
xinit = [0.890; 1.096; 0.159; 426; 0; 
    0.928; 1.142; 0.159; 426; 0;
    0.954];

[Ts, xsize, ~, usize, ysize, uoff1, uoff2, ud] = const_sim();
u = [uoff1; uoff2; 0] + [Td(1); 0; 0; u_rec(1);Td(2); 0; 0; u_rec(2); ud];

% initial derivatives
[f1,m1] = get_comp_deriv(xinit(1:xsize),[u(1:usize); Pin; P_D],1);
[f2,m2]= get_comp_deriv(xinit(xsize+1:2*xsize),[u(usize+1:2*usize); P_D; Pout],1);
ftank = get_tank_deriv(P_D,[m1+m2; ud; xinit(xsize+1)]);

% linearize system
[Ac,Bc,Ccorig] = linearize_tank(xinit, u);
Cc = [Ccorig([2,4],:); Ccorig(1,:)-Ccorig(3,:); Ccorig(5,:)];

% disturbances
delta_x = ones(N,1)*[repmat([0.01, 0.01, 0.05, 10, 0.05],1,2),0.02];
delta_x = (rand(2*xsize+1,N)-0.5).*delta_x';

delta_u = 0*ones(N,1)*[repmat([0.01, 0, 0, 0.01],1,2),0];
delta_u = (rand(2*usize+1,N)-0.5).*delta_u';

f = zeros(N,2*xsize+1);
fest = f;

for i=1:N
    P_D = P_D_init + delta_x(end,i);
    [fe1,m1] = get_comp_deriv(xinit(1:xsize)+delta_x(1:xsize,i),[u(1:usize)+delta_u(1:usize,i); Pin; P_D],1);
    [fe2,m2] = get_comp_deriv(xinit(xsize+1:2*xsize)+delta_x(xsize+1:2*xsize,i),[u(usize+1:2*usize)+delta_u(usize+1:2*usize,i); P_D; Pout],1);
    fetank = get_tank_deriv(xinit(end)+delta_x(end,i), [m1+m2; u(end)+delta_u(end,i); xinit(xsize+1)+delta_x(xsize+1,i)]);
    
    f(i,:) = [fe1-f1; fe2-f2; fetank-ftank];
    
    fest(i,:) = Ac*delta_x(:,i) + Bc*delta_u([1,4,usize+1,usize+4],i);
    
end
close all

for i=1:11
    figure; 
    plot(1:N,[f(:,i),fest(:,i)])
end


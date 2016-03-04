% N = 200;
% rec_open = zeros(N,1);
% % Td = 0.03*ones(N,1);
% Td = linspace(0,0.05,N)';
% % Td = rec_open;
% Inflow_opening = 0.405;
% Outflow_opening = 0.393;
% 
% Ts = 0.05;
% 
% t = 0:Ts:N*Ts;
% 
% x_init_lin = [0.898
%       1.126
%       0.15
%       439.5
%       0];
% u_init = [0.304+Td(1),Inflow_opening,Outflow_opening,rec_open(1)]';
% 
% 
% x = zeros(N+1,5);
% dx = x;

% x(1,:) = x_init_lin';
% x(1,2) = x(1,2)-1.126;
% f = get_comp_deriv(x_init_lin,u_init);
% x(2,:) = x_init_lin'+f*Ts;
% x(1,:) = x_out.signals.values(1,:);
% x(1,2) = x(1,2) + 1.126;
% 
% f = get_comp_deriv(x(1,:)',[0.304+u(1,1),Inflow_opening,Outflow_opening,u(1,2)]);
% [Ac,Bc,Cc] = get_linearized_matrices(x(1,:)',u(1,:)');
% 
% x(2,:) = x(1,:) + ((Ts + Ts/2*Ac + Ts^2/6*Ac^2 + Ts^3/24*Ac^3)*f')';
% x(1:2,:) = x_out.signals.values(1:2,:);
% 
% x(1,:) = x_out.signals.values(1,:);
% f = get_comp_deriv(x(1,:)',u_init);
% [Ac,Bc,Cc] = get_linearized_matrices(x(1,:)',u(1,:)');
% [A,B,C] = discretize_rk4(Ac,Bc,Cc,Ts);
% x(2,:) = x(1,:)+f*Ts;
% y(1,1) = p_out.signals.values(1,1);
% y(1,2) = SD(1);
% 
% y(2,:) = y(1,:) + (C*(x(2,:)-x(1,:))')';
% 
% u = [Td,rec_open];
% 
% for i=2:N
% %     x(i,2) = x(i,2)-1.126;
%     [Ac,Bc,Cc] = get_linearized_matrices(x(i,:)',u(i-1,:)');
%     [A,B,C] = discretize_rk4(Ac,Bc,Cc,Ts);
%     f = get_comp_deriv(x(i,:)',[0.304+u(i-1,1),Inflow_opening,Outflow_opening,u(i-1,2)]);
%     f2(i,:) = f;
%     dx(i-1,:) = ((Ts*eye(size(Ac)) + Ts^2/2*Ac + Ts^3/6*Ac^2 + Ts^4/24*Ac^3)*f')';
%     x(i+1,:) = (B*(u(i,:)-u(i-1,:))')' + x(i,:) + dx(i-1,:);
%     y(i+1,:) = (C*(x(i+1,:)-x(i,:))')'+y(i,:);
% %     x(i+1,:) = x(i,:) + dx(i-1,:);
% %     x(i+1,:) = x(i,:) + (x(i,:)-x(i-1,:));
% end

%%
clear all
N = 200;
rec_open = zeros(N,1);
rec_open = 0.5 *ones(N,1);
% Td = 0.03*ones(N,1);
% Td = [linspace(0,0.01,N/2),linspace(0.01,0,N/2)]';
Td = 0.1*ones(N,1);
Inflow_opening = 0.405;
Outflow_opening = 0.393;

Ts = 0.05;

x_init_lin = [0.898
      1.126
      0.15
      439.5
      0];
%   load config.mat
%   
% rec_open = u_apply(2,:)';
% Td = u_apply(1,:)';
% N = length(rec_open);
% x_init_lin = x_init;
  
% x_init_lin = [0.899; 1.125; 0.1512; 440.7; 0];

u_init = [0.304+Td(1),Inflow_opening,Outflow_opening,rec_open(1)]';

if ~(exist('p_out','var'))
    sim('compressor.mdl')
end
% t = p_out.time;
t = 0:Ts:(N-1)*Ts;
y = [p_out.signals.values,SD];
yr = interp1(p_out.time,y,t)';

u = [Td,rec_open]';

n_delay = [0;20];

% [Ac,Bc,Cc] = get_linearized_matrices(x_init_lin,u(:,1));
% [A,B,C] = discretize_rk4(Ac,Bc,Cc,Ts);
load sys.mat

[Aaug,Baug,Caug] = get_augmented_matrices(A,B,C,n_delay);

orig_xsize = 5;
dsize = 2;
ysize = 2;
usize = 2;
xsize = length(Aaug);

x = zeros(5,length(t));
dx = zeros(xsize,length(t));
xhat = dx;
x(1:5,1) = x_init_lin;
x(1:5,2) = x_init_lin;

% expectation of output disturbance
Qn = eye(orig_xsize);
Rn = eye(dsize);

% state disturbance to state matrix
G = [zeros(xsize-dsize,orig_xsize);
    zeros(dsize,1), eye(dsize)*10, zeros(dsize,orig_xsize-dsize-1)];

Hkalman = [zeros(dsize,orig_xsize-dsize), eye(dsize)];

sys_kalman = ss(Aaug,[Baug G], Caug, [zeros(ysize,usize), Hkalman], Ts);

Caug_init = Caug;

[KEST,L,P,M,Z] = kalman(sys_kalman,Qn,Rn); % get observer values

yt(:,1) = yr(:,1);

xaug = dx;
dxaug = dx;
% udele = zeros(xsize-5,length(t));
xaug = zeros(5,length(t));
xaug(1:5,1) = x_init_lin;
xaug(1:5,2) = x_init_lin;

for i=2:length(t)
    dxaug(:,i) = dxaug(:,i) + M*(yr(:,i)-yr(:,i-1) - Caug*(dxaug(:,i)));
    xaug(:,i) = xaug(:,i-1) + dxaug(1:5,i);
%     x(:,i) = x(:,i) + M*(yr(:,i) - Caug*(x(:,i)-x(:,i-1)));
%     dx(:,i) = dx(:,i) + M*(yr(:,i) - yr(:,i-1) - Caug*dx(:,i));
%     x(:,i) = x(:,i-1) + dx(1:5,i);
    
%     yt(:,i) = Caug*(x(:,i)-x(:,i-1))+yt(:,i-1);
%     ydiff(:,i) = Caug*(x(:,i)-x(:,i-1))-(yr(:,i)-yr(:,i-1));
%     
%     [Ac,Bc,Cc] = get_linearized_matrices(xaug(1:5,i),[u(1,i-1); dxaug(6,i)]);
%     f = get_comp_deriv(xaug(1:5,i),[0.304+u(1,i-1),Inflow_opening,Outflow_opening,dxaug(6,i)])';

    [Ac,Bc,Cc] = get_linearized_matrices(xaug(1:5,i),[u(1,i-1); u(2,i-1)]);
    f = get_comp_deriv(xaug(1:5,i),[0.304+u(1,i-1),Inflow_opening,Outflow_opening,u(2,i-1)])';
    [A,B,C,fd] = discretize_rk4(Ac,Bc,Cc,f,Ts);
    
    [Aaug,Baug,Caug] = get_augmented_matrices(A,B,C,n_delay);
    
    xlin = dxaug(:,i);
    xlin(6) = u(2,i-1);
    xlin(7:end) = 0;
    
    dxaug(:,i+1) = Baug*([u(1,i)-u(1,i-1); u(2,i)]) + Aaug*(dxaug(:,i)-xlin) + [fd; zeros(size(dxaug,1)-5,1)];
%     dxaug(1:5,i) - fd
%     xaug(1:5,i+1) = xaug(1:5,i+1)+xaug(1:5,i)+fd;
%     dxaug(1:5,i+1) = dxaug(1:5,i+1)+fd;%+Aaug(1:5,6)*xaug(6,i+1);
%     dxaug(6:end,i+1) = dxaug(6:end,i+1)+Aaug(6:end,6:end)*dxaug(6:end,i);
    
    
    


end

% 
% xr = x_out.signals.values;
% for i=2:length(x_out.time)
%     fr(i,:) = get_comp_deriv(xr(i,:)',[0.304+u(min(i-1,N),1),Inflow_opening,Outflow_opening,u(min(i-1,N),2)]);
%     fr2(i,:) = (xr(i,:)-xr(i-1,:))/(x_out.time(i)-x_out.time(i-1));
% end




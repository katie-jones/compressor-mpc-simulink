x_init_lin = [0.899; 1.125; 0.1512; 440.7; 0];

x_init_lin = [0.898
      1.126
      0.15
      439.5
      0];
  
  
upast = [0.1;0];

udiff = [0;0.1];

[orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m] = constants();

xinit = [x_init_lin; upast(2)*ones(n_delay(2),1); zeros(dsize,1)];

[A,B,C,H,Ga,Gb,Gc,dx,Sx,Su,Sf] = get_qp_matrices(x_init_lin,upast);

dx = zeros(size(dx));

x0 = xinit - [x_init_lin; upast(2)*ones(n_delay(2),1); zeros(2,1)];

y = Sx*x0 + Su * repmat(udiff,m,1) + Sf * dx;

y = reshape(y,2,p);

y2 = zeros(size(y));

x = x0;
for i=1:p
    x = A*x + B*udiff + dx;
    y2(:,i) = C*x;
end

%% test QP
clear all
addpath('../call_qpoases_m')
addpath('../call_qpoases_m/qpoases3/interfaces/matlab/')
x_init_lin = [0.898; 1.126; 0.15; 439.5; 0];
% x_init_lin = [0.899; 1.125; 0.1512; 440.7; 0];

[C,~,UB,LB,b,Ain] = get_constant_matrices();
[orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m] = constants();

if ~(exist('upast','var'))
    upast = [0; 0];
    xinit = [x_init_lin; upast(2)*ones(n_delay(2),1); zeros(2,1)];
    dxaug = zeros(xsize,1);
end
% yref = [1.126; 4.3];
% y = [1.126;4.3];
yref = [0;0];
y = [0;0];

N = 100;
t = 0:Ts:(N-1)*Ts;

for k=1:N

xinit(1:5) = xinit(1:5) + dxaug(1:5);
xinit(6:end) = dxaug(6:end);

x_out(:,k) = xinit;

y = y + C*(dxaug);

y_out(:,k) = y;


% Get reference and input vector
Yref = zeros(p*ysize,1);
Upast = zeros(m*usize,1);
for i=1:p
    Yref(1+(i-1)*ysize:i*ysize,1) = yref - y;
end
for i=1:m
    Upast(1+(i-1)*usize,1) = upast(1);
%     Upast(i*usize,1) = upast(2);
    Upast(i*usize,1) = xinit(6);
end

% Update bounds
LBa = LB - Upast;
UBa = UB - Upast;

[A,B,C,H,Ga,Gb,Gc,df0,Sx,Su,Sf,UWT] = get_qp_matrices(xinit,[upast(1); xinit(6)]);

% Calculate initial state
dx0 = [zeros(6,1); xinit(7:end-2)-xinit(6); xinit(end-1:end)];

% Gradient vector
f = df0'*Gc - Yref'*Ga + dx0'*Gb + Upast'*UWT;

% Solve QP
lbA = -inf(size(Ain,1),1);
[usol,cost,exitflag,it,lambda] = qpOASES(H,f',Ain,LBa,UBa,lbA,b);

if (exitflag~=0)
    usol=quadprog(H,f,[eye(4); -eye(4); Ain],[UBa; -LBa; b]);
    k
end
exitflag_out = exitflag(1,1);
it_out = it(1,1);

% Apply next input
delta_u = usol(1:usize,1);
% u_new = [upast(1); dxaug(6)] + delta_u;
u_new(:,k) = [upast(1); dxaug(6)] +delta_u;

% value of x about which we linearized
xlin = dxaug;
% xlin(6) = upast(2);
xlin(7:end) = 0;

dxaug = B*[delta_u(1); u_new(2,k)] + A*(dxaug-xlin) + df0;
upast = u_new(:,k);

if k==50
%     keyboard
end

end

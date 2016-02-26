function [u_new,exitflag_out,it_out,x_out]  = controller_lin(u)
%#eml
% persistent C_new M LB UB Ga Gb Gc H Ain b A B upast xinit
persistent M LB UB Ain b upast xinit dxaug y_old C %x_init_lin
eml.extrinsic('disp');
eml.extrinsic('toc');
eml.extrinsic('tic');
eml.extrinsic('num2str');
eml.extrinsic('keyboard');

p = 100;
m = 2;

if isempty(M)
    tic
    option = 1;
%     [H,C_new,M,UB,LB,b,Ain,Ga,Gb,Gc,xinit,upast] = get_matrices(option);
    [C,M,UB,LB,b,Ain,xinit,upast] = get_constant_matrices();
    dxaug = zeros(size(xinit));
    y_old = [0; 0];
    time = toc;
    disp(['Preparation time: ' num2str(time)]);
end

y = u(1:2);
y(1) = y(1) + 1.126;
y(2) = y(2) + 4.3;

% correct for normalization in compressor model
yref = u(3:4);
yref(1) = yref(1) + 1.126;
yref(2) = yref(2) + 4.3;

% Estimate state info
% dxaug = dxaug + M*(y - y_old - C*(dxaug));
xinit(1:5) = xinit(1:5) + dxaug(1:5);
xinit(6:end) = dxaug(6:end);
x_out = xinit;

% Get matrices - linearize about delayed rec.open input
% [A,B,C,H,Ga,Gb,Gc,dx] = get_qp_matrices(xinit,[upast(1); dxaug(6)]);
[A,B,C,H,Ga,Gb,Gc,df0] = get_qp_matrices(xinit,upast);

% Get reference and input vector
ysize = 2;
usize = 2;
Yref = zeros(p*ysize,1);
Upast = zeros(m*usize,1);
for i=1:p
    Yref(1+(i-1)*ysize:i*ysize,1) = yref - y;
end
for i=1:m
    % use delayed recycle opening as reference input
    Upast(1+(i-1)*usize,1) = upast(1);
    Upast(i*usize,1) = upast(2);
%     Upast(i*usize,1) = xinit(6);
end

% Update bounds
LBa = LB - Upast;
UBa = UB - Upast;

% Calculate initial state
dx0 = [zeros(5,1); xinit(6:end-2)-upast(2); xinit(end-1:end)];

% Gradient vector
f = df0'*Gc - Yref'*Ga + dx0'*Gb;

% Solve QP
lbA = -inf(size(Ain,1),1);
[usol,cost,exitflag,it,lambda] = call_qpoases_m(H,f,Ain,LBa,UBa,lbA,b);
exitflag_out = exitflag(1,1);
it_out = it(1,1);

% Apply next input
delta_u = usol(1:usize,1);
% u_new = [upast(1); dxaug(6)] + delta_u;
u_new = upast + delta_u;

% value of x about which we linearized
xlin = dxaug;
xlin(6) = upast(2);
xlin(7:end) = 0;

dxaug = B*[u_new(1)-upast(1);u_new(2)] + A*(dxaug-xlin) + df0;

upast = u_new;


% % p=200;
% % m=8;
% p = 100;
% m = 2;
% 
% if isempty(M)
%     tic
%     option = 1;
% %     [H,C_new,M,UB,LB,b,Ain,Ga,Gb,Gc,xinit,upast] = get_matrices(option);
%     [C,M,UB,LB,b,Ain,xinit,upast,deltax] = get_constant_matrices();
% %     x_init_lin = [0.8980 1.1260 0.1500 439.5000 0]';
%     time = toc;
%     disp(['Preparation time: ' num2str(time)]);
% end
% 
% % tic
% 
% y=u(1:2);
% yref=u(3:4);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % estimate the state information
% 
% 
% % xinit  = xinit + M*(y - C_new*xinit);
% deltax = deltax + M*(y - C*deltax);
% xinit = deltax + xinit;
% x_out = xinit;
% 
% [A,B,C,H,Ga,Gb,Gc] = get_qp_matrices(xinit,upast);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% % upper and lower bounds wrt past input
% Upast = zeros(m*2,1);
% 
% for i=1:m
%     Upast(1+(i-1)*2:i*2,1)=upast;
% end
% 
% LBa=LB-Upast;
% UBa=UB-Upast;
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % generate the output reference vector
% 
% R = zeros(p*2,1);
% 
% for i=1:p
% %     R(1+(i-1)*2:i*2,1)=yref - y;
%     R(1+(i-1)*2:i*2,1)=y - yref;
% end
% 
% % generate the gradient vector
% 
% % f = deltax'*Gb+Upast'*Gc-R'*Ga;
% f = deltax'*Gb-R'*Ga;
% 
% % xinit1 = xinit;
% % upast1 = upast;
% 
% lbA = -inf(size(Ain,1),1);
% ubA = b;
% [usol,cost,exitflag,it,lambda] = call_qpoases_m(H,f,Ain,LBa,UBa,lbA,ubA);
% exitflag_out = exitflag(1,1);
% it_out = it(1,1);
% dunext=usol(1:2,1);
% 
% deltax= A*deltax + B*(dunext);
% 
% u_new = dunext + upast;
% % u_new = zeros(2,1);
% 
% upast = u_new;

% time2 = toc;
% disp(['Normal Execution Time: ' num2str(time2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% solve the QP

% usol=quadprog(H_qp,f,Ain,b,[],[],LBa,UBa);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % extract the control moves for next time step
% 
% 
% dunext=usol(1:2,1);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % initial conditions are updated
% 
% 
% xinit= A*xinit + B*(dunext+upast);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % control moves are sent to the plant
% 
% upast=dunext+upast;
% 
% sys=upast;

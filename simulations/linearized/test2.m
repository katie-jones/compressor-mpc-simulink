clear all
x_init_lin = [0.899; 1.125; 0.1512; 440.7; 0];
x_init_lin2 = [
    0.8990
    1.1250
    0.1512
  440.6442
         0];

% x_init_lin = [     0.9293    1.0702    0.1633  357.3388    0.0403]';

[C,~,UB,LB,b,Ain] = get_constant_matrices();
[orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m] = constants();

u1 = [-0.008; 0.0003];
u2 = [-0.0439; 0.0251];

u_apply = [u1,repmat(u2,1,p-1)];

t = 0:0.05:(p-1)*0.05;

rec_open = u_apply(2,:)';
Td = u_apply(1,:)';

yref = [0;0];
% y = [0;0];
y = C(:,1:5)*(x_init_lin - x_init_lin2);


xinit = [x_init_lin; zeros(n_delay(2)+dsize,1)];

[A,B,C,H,Ga,Gb,Gc,df0,Sx,Su,Sf,UWT] = get_qp_matrices(xinit,u_apply(:,1));
dx0 = [zeros(5,1); xinit(6:end-2)-u_apply(2,1); xinit(end-1:end)];

y2 = Su*reshape(u_apply(:,1:m),m*usize,1) + Sx*dx0 + Sf*df0;
y2 = reshape(y2,2,p);
y2 = y2 + repmat(y,1,p);

sim('compressor')

figure; 
plot(p_out.time,SD*100-4.3,t,y2(2,:))
xlim([0 p*0.05])
function [usol,flag,it] = run_qp(xinit,upast,yref,y)

[C,~,UB,LB,b,Ain] = get_constant_matrices();

[orig_xsize,ysize,dsize,usize,n_delay,xsize,Ts,p,m] = constants();

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


end
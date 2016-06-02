function [A,B,C,H,Ga,Gb,Gc,dx,Sx,Su,Sf,UWT] = get_qp_matrices(xinit,upast)
%#eml
%% Constants
Inflow_opening = 0.405;
Outflow_opening = 0.393;

[xsize,ysize,dsize,usize,n_delay,xsize_full,Ts,p,m,UWT,YWT] = const_mpc();

%% Linearized system
u = [0.304+upast(1),Inflow_opening,Outflow_opening,upast(2),0]; % include inflow/outflow openings
[Ac,Bc,Cc] = get_linearized_matrices(xinit,u);

f = get_comp_deriv(xinit(1:5),u,1);

[Ainit,Binit,Cinit,dx2] = discretize_rk4(Ac,Bc,Cc,f,Ts);

%% Define augmented system
% delay of n_delay time steps in 2nd component of u
% constant output disturbance

Adelay2 = [zeros(n_delay(2)-1,1), eye(n_delay(2)-1); zeros(1,n_delay(2))]; % delay component of A
Adist = eye(dsize); % disturbance component of A

if (n_delay(1)==0)
    A = [Ainit, Binit(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(2),xsize), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(2)), Adist];

    B = [Binit(:,1), zeros(xsize,1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1];
        zeros(dsize,2)];
else
    Adelay1 = [zeros(n_delay(1)-1,1), eye(n_delay(1)-1); zeros(1,n_delay(1))]; % delay component of A

    A = [Ainit, Binit(:,1), zeros(xsize,n_delay(1)-1), Binit(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(1),xsize), Adelay1, zeros(n_delay(1),n_delay(2)), zeros(n_delay(1),dsize);
        zeros(n_delay(2),xsize), zeros(n_delay(2),n_delay(1)), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(1)), zeros(dsize,n_delay(2)), Adist];

    B = [zeros(xsize,usize);
        [zeros(n_delay(1)-1,1); 1], zeros(n_delay(1),1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1); 1];
        zeros(dsize,usize)];
end

% disturbance to output matrix
Cd = eye(ysize,dsize); 

C = [Cinit, zeros(ysize,sum(n_delay)), Cd];

%% Define correct global variables
xsize = xsize_full;

dx = [dx2; zeros(xsize-5,1)];


%% Define system matrices

% Y = Su*U + Sx*X

% Pre-compute multiples of C*A^(i-1)
CxA = zeros(ysize,xsize,p+1);
CxA(:,:,1) = C;
for i=2:p+1
    CxA(:,:,i) = CxA(:,:,i-1)*A;
end

Su = zeros(ysize*p,usize*m);


for i=1:p
    for j=1:i
        % for first m inputs, make new columns
        if j<=m
            Su(1+(i-1)*ysize:i*ysize,1+(j-1)*usize:j*usize) = CxA(:,:,i-j+1)*B;
            
        % m+1:p inputs are the same as input m
        else
            toadd = CxA(:,:,i-j+1)*B;
            for k=1:ysize
                Su(k+(i-1)*ysize,1+(m-1)*usize:m*usize) = Su(k+(i-1)*ysize,1+(m-1)*usize:m*usize) + toadd(k,:);
            end
        end
    end
end


Sx = zeros(ysize*p,xsize);
Sf = Sx;
for i=1:p
    for j=1:ysize
        Sx(j+(i-1)*ysize,:) = CxA(j,:,i+1);
    end
end


Sf(1:ysize,:) = C;
for i=2:p
    for j=1:ysize
        Sf(j+(i-1)*ysize,:) = Sf(j+(i-2)*ysize,:) + CxA(j,:,i);
    end
end


%% Calculate QP matrices

H = Su'*YWT*Su + UWT;

Ga = YWT*Su;
Gb = Sx'*YWT*Su;
Gc = Sf'*YWT*Su;



end



% Discretize system given by A,B,C using 4th order runge kutta
function [Ad,Bd,Cd,fd] = discretize_rk4(A,B,C,f,Ts)
%#eml
xsize = 5;

A2 = A*A;
A3 = A2*A;

Acom = Ts*eye(xsize) + Ts^2/2*A + Ts^3/6*A2 + Ts^4/24*A3;

Ad = eye(xsize) + Acom*A;
Bd = Acom*B;
Cd = C;
fd = Acom*f(:);

end

function [A,B,C,H,Ga,Gb,Gc,dx,Sx,Su,Sf,UWT] = get_qp_matrices(xinit,upast)

[xsize_comp, xsize, ~, ~, ~, ysize, uoff1, uoff2, ud] = const_sim();
[Ts,n_delay,dsize,usize,p,m,UWT,YWT] = const_mpc2();

x1 = xinit(1:xsize_comp);
x2 = xinit(xsize_comp+1:2*xsize_comp);
pd = xinit(2*xsize_comp+1);

u1 = uoff1 + [upast(1); 0; 0; upast(2); 0];
u2 = uoff2 + [upast(usize+1); 0; 0; upast(usize+2); 0];

u1(end) = pd; % give output pressure as last input
u2(end) = pd;

u = [u1; u2; ud];

[Ac,Bc,Cc] = linearize_tank(xinit,u);

f1 = get_comp_deriv(x1,u1,0);
f2 = get_comp_deriv(x2,u2,0);
ftank = get_tank_deriv(xinit,u);

[Ainit,Binit,Cinit,dx2] = discretize_rk4(Ac,Bc,Cc,[f1; f2; ftank],Ts);

%% Augment matrices

Adelay2 = [zeros(n_delay(2)-1,1), eye(n_delay(2)-1); zeros(1,n_delay(2))]; % delay component of A
Adist = eye(2*dsize); % disturbance component of A


A = [Ainit, Binit(:,2), zeros(xsize,n_delay(2)-1), Binit(:,usize+2), zeros(xsize,n_delay(2)-1), zeros(xsize,2*dsize);
    zeros(n_delay(2),xsize), Adelay2, zeros(n_delay(2)), zeros(n_delay(2),2*dsize);
    zeros(n_delay(2),xsize), zeros(n_delay(2)), Adelay2, zeros(n_delay(2),2*dsize);
    zeros(2*dsize,xsize), zeros(2*dsize,2*n_delay(2)), Adist];

B = [Binit(:,1), zeros(xsize,1), Binit(:,usize+1), zeros(xsize,1);
    zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1], zeros(n_delay(2),1), zeros(n_delay(2),1);
    zeros(n_delay(2),1), zeros(n_delay(2),1), zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1];
    zeros(2*dsize,2*usize)];
    
Cdist = eye(ysize,2*dsize);
C = [Cinit(1:ysize,:), zeros(ysize,2*sum(n_delay)), Cdist];

dx = [dx2; zeros(2*sum(n_delay)+2*dsize,1)];



%% Define system matrices

% Y = Su*U + Sx*X
xsize = xsize + 2*sum(n_delay) + 2*dsize;

% Pre-compute multiples of C*A^(i-1)
CxA = zeros(ysize,xsize,p+1);
CxA(:,:,1) = C;
for i=2:p+1
    CxA(:,:,i) = CxA(:,:,i-1)*A;
end

Su = zeros(ysize*p,2*usize*m);


for i=1:p
    for j=1:i
        % for first m inputs, make new columns
        if j<=m
            Su(1+(i-1)*ysize:i*ysize,1+(j-1)*2*usize:j*2*usize) = CxA(:,:,i-j+1)*B;
            
        % m+1:p inputs are the same as input m
        else
            toadd = CxA(:,:,i-j+1)*B;
            for k=1:ysize
                Su(k+(i-1)*ysize,1+(m-1)*2*usize:m*2*usize) = Su(k+(i-1)*ysize,1+(m-1)*2*usize:m*2*usize) + toadd(k,:);
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

A2 = A*A;
A3 = A2*A;

Acom = Ts*eye(size(A)) + Ts^2/2*A + Ts^3/6*A2 + Ts^4/24*A3;

Ad = eye(size(A)) + Acom*A;
Bd = Acom*B;
Cd = C;
fd = Acom*f(:);

end


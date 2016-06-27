function [A,B,C,dx,H1,H2,f0_1,f0_2,Gd1,Gd2] = get_qp_matrices(xinit,upast,dyref,UWT,YWT,~)

%% Constants
[Ts,xsize_comp, xsize, ~, ysize, uoff1, uoff2] = const_sim();
[n_delay,dsize,usize,p,m] = const_mpc();
[~,Pin,Pout] = const_flow();

x1 = xinit(1:xsize_comp);
x2 = xinit(xsize_comp+1:2*xsize_comp);

u1 = [uoff1 + [upast(1); 0; 0; upast(2)]; Pin; x2(1)];
u2 = [uoff2 + [upast(usize+1); 0; 0; upast(usize+2)]; -1; Pout;];

[f1,m_out1] = get_comp_deriv(x1,u1,1);
[f2,~] = get_comp_deriv(x2,[u2; m_out1],1);

u = [u1; u2];

%% Linearization

[Ac,Bc,Ccorig] = linearize_serial(xinit, u);
Cc = Ccorig(1:ysize,:);

[Ainit,Binit,Cinit,dx2] = discretize_rk4(Ac,Bc,Cc,[f1; f2],Ts);

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

% derivative at linearization point
dx = [dx2; zeros(2*sum(n_delay)+2*dsize,1)];



%% Define system matrices

% Y = Su*U + Sx*X
xtotalsize = xsize + 2*sum(n_delay) + 2*dsize;

% Pre-compute multiples of C*A^(i-1)
CxA = zeros(ysize,xtotalsize,p+1);
CxA(:,:,1) = C;
for i=2:p+1
    CxA(:,:,i) = CxA(:,:,i-1)*A;
end

Su1 = zeros(ysize*p,usize*m);
Su2 = Su1;

B1 = B(:,1:usize);
B2 = B(:,usize+1:end);

for i=1:p
    for j=1:i
        % for first m inputs, make new columns
        if j<=m
            Su1(1+(i-1)*ysize:i*ysize,1+(j-1)*usize:j*usize) = CxA(:,:,i-j+1)*B1;
            Su2(1+(i-1)*ysize:i*ysize,1+(j-1)*usize:j*usize) = CxA(:,:,i-j+1)*B2;
            
        % m+1:p inputs are the same as input m
        else
            toadd1 = CxA(:,:,i-j+1)*B1;
            toadd2 = CxA(:,:,i-j+1)*B2;
            for k=1:ysize
                Su1(k+(i-1)*ysize,1+(m-1)*usize:m*usize) = Su1(k+(i-1)*ysize,1+(m-1)*usize:m*usize) + toadd1(k,:);
                Su2(k+(i-1)*ysize,1+(m-1)*usize:m*usize) = Su2(k+(i-1)*ysize,1+(m-1)*usize:m*usize) + toadd2(k,:);
            end
        end
    end
end


Sx = zeros(ysize*p,xtotalsize);
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

%% Calculate QP matrices for two compressors
% deltax0 (augmented)
deltax0 = [zeros(xsize,1); xinit(xsize+1:xsize+n_delay(2))-upast(2); xinit(xsize+n_delay(2)+1:xsize+2*n_delay(2))-upast(usize+2); xinit(end-2*dsize+1:end)];

% reference vector
Yref = zeros(p*ysize,1);
for i=1:p
    Yref(1+(i-1)*ysize:i*ysize,1) = dyref;
end

% Quadratic term for each compressor
H1 = Su1'*YWT*Su1 + UWT;
H2 = Su2'*YWT*Su2 + UWT;

% Cross terms
Ga1 = YWT*Su1; % yref cross term
Gb1 = Sx'*YWT*Su1; % deltax0 cross term
Gc1 = Sf'*YWT*Su1; % dx (derivative at lin. pt.) cross term
Gd1 = Su2'*YWT*Su1; % u_other cross term

Ga2 = YWT*Su2;
Gb2 = Sx'*YWT*Su2;
Gc2 = Sf'*YWT*Su2;
Gd2 = Su1'*YWT*Su2;

% Gradient vector
f0_1 = dx'*Gc1 - Yref'*Ga1 + deltax0'*Gb1;
f0_2 = dx'*Gc2 - Yref'*Ga2 + deltax0'*Gb2;


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


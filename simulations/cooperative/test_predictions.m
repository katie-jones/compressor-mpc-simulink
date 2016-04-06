% Demonstrate effect of Td on surge distances for various conditions
function test_predictions()
addpath('../common')
addpath('../parallel_common/')

save_plots=0;

%% Constants
[Ts, ~, ~, ~, ysize] = const_sim();
[~,~,ucontrolsize,p,m] = const_mpc();

P_D = 1.12;

x_init_lin = [0.916; 1.145; 0.152; 440; 0];

xinit = [x_init_lin; x_init_lin; P_D];

uinit = zeros(2*ucontrolsize,1);

%% Predictions
Su1 = get_prediction_matrix(xinit,uinit);

du = repmat([0.01; 0],m,1); % step on torque, no recycle valve

y1 = reshape(Su1*du,ysize,p)'; % prediction of compressor 1

t = 1:p;

SD_sum = zeros(p,2);

SD_sum(1,:) = [0 0];
for i=2:p
    SD_sum(i,:) = SD_sum(i-1,:) + Ts*y1(i,1:2);
end
%%

fig=figure; 
set(fig,'units','centimeters','position',[35 5 28 15])
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',2)
subplot(1,2,1)
plot(t,y1(:,1:2));
grid on
title('Predicted surge distance')
xlabel({'Horizon length (m)'; 'x = [0.92 1.15 0.15 440 0]'''})
ylabel('Relative surge distance [%]')
legend('SD_1','SD_2','location','southwest')

subplot(1,2,2)
plot(t,SD_sum.^2);

ylim([0 0.25])
xlim([0 p])
grid on
title('Cumulative effect of predicted surge distance')
xlabel({'Horizon length (m)'; '\Delta T_{d,1}=0.01, p=2'})
ylabel('Cumulative effect of surge distance')

if save_plots
    fig=printplot(fig);
    saveas(fig,'Td_effect.pdf')
end

set(0,'defaultlinelinewidth',1)
end


function Su1 = get_prediction_matrix(xinit,upast)

%% Constants
[Ts,xsize_comp, xsize, ~, ysize, uoff1, uoff2, ud] = const_sim();
[n_delay,dsize,usize,p,m] = const_mpc();

x1 = xinit(1:xsize_comp);
x2 = xinit(xsize_comp+1:2*xsize_comp);
pd = xinit(2*xsize_comp+1);

u1 = uoff1 + [upast(1); 0; 0; upast(2); 0];
u2 = uoff2 + [upast(usize+1); 0; 0; upast(usize+2); 0];

u1(end) = pd; % give output pressure as last input
u2(end) = pd;

u = [u1; u2; ud];

%% Linearization

[Ac,Bc,Cc] = linearize_tank(xinit,u);

f1 = get_comp_deriv(x1,u1,0);
f2 = get_comp_deriv(x2,u2,0);
ftank = get_tank_deriv(xinit,u);

[Ainit,Binit,Cinit] = discretize_rk4(Ac,Bc,Cc,[f1; f2; ftank],Ts);

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


%% Define system matrices

% Y = Su*U + Sx*X
xsize = xsize + 2*sum(n_delay) + 2*dsize;

% Pre-compute multiples of C*A^(i-1)
CxA = zeros(ysize,xsize,p+1);
CxA(:,:,1) = C;
for i=2:p+1
    CxA(:,:,i) = CxA(:,:,i-1)*A;
end

Su1 = zeros(ysize*p,usize*m);

B1 = B(:,1:usize);

for i=1:p
    for j=1:i
        % for first m inputs, make new columns
        if j<=m
            Su1(1+(i-1)*ysize:i*ysize,1+(j-1)*usize:j*usize) = CxA(:,:,i-j+1)*B1;
            
        % m+1:p inputs are the same as input m
        else
            toadd1 = CxA(:,:,i-j+1)*B1;
            for k=1:ysize
                Su1(k+(i-1)*ysize,1+(m-1)*usize:m*usize) = Su1(k+(i-1)*ysize,1+(m-1)*usize:m*usize) + toadd1(k,:);
            end
        end
    end
end

end


% Discretize system given by A,B,C using 4th order runge kutta
function [Ad,Bd,Cd,fd] = discretize_rk4(A,B,C,f,Ts)

A2 = A*A;
A3 = A2*A;

Acom = Ts*eye(size(A)) + Ts^2/2*A + Ts^3/6*A2 + Ts^4/24*A3;

Ad = eye(size(A)) + Acom*A;
Bd = Acom*B;
Cd = C;
fd = Acom*f(:);

end

function [n_delay,dsize,usize,p,m] = const_mpc()
n_delay = [0; 40];
dsize = 2;
usize = 2;
p = 100;
m = 2;

end

% Function to linearize discharge tank about point xbar, ubar

function [A, B, C] = linearize_serial(x,u)
%#eml

%% Constants

xsize = 5;
usize = 6;
ysize = 2;

x1 = x(1:xsize);
x2 = x(xsize+1:2*xsize);

u1 = u(1:usize);
u2 = u(usize+1:2*usize);

uoutput1 = u1(3);

[~,~,~,~,~,~,D2] = comp_coeffs();

[SpeedSound,~,~,V1,V2] = const_flow();

%% Compressor dynamics
[A1,B1,C1] = get_linearized_matrices(x1,u1);

[A2,B2,C2] = get_linearized_matrices(x2,u2);

%% Fix dp12/p12
A2(1,1) = -getValveDerivative(x1(2), x2(1), SpeedSound, V1, uoutput1, D2);

%% Cross terms
A12 = zeros(xsize); % effect of compressor 1 on compressor 2
A21 = A12; % effect of compressor 2 on compressor 1
A12(1,2) = getValveDerivative(x1(2), x2(1), SpeedSound, V1, uoutput1, D2);
A21(2,1) = getValveDerivative(x1(2), x2(1), SpeedSound, V2, uoutput1, D2);

%% Output system
A = [A1, A21;
    A12, A2];

B = [B1, zeros(xsize,2);
    zeros(xsize,2), B2];

% Outputs: p21,SD1,p22,SD2,PD
% C = [C1, zeros(ysize,xsize), zeros(ysize,1);
%     zeros(ysize,xsize), C2, zeros(ysize,1);
%     zeros(1,2*xsize), 1];

% Outputs: p21,SD1,p22,SD2,mc1-mc2
C = [C1, zeros(ysize,xsize);
    zeros(ysize,xsize), C2;
    [0,0,1,0,0],[0,0,-1,0,0]];



end



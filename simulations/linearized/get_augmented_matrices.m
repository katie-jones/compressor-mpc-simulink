function [Aaug,Baug,Caug] = get_augmented_matrices(A,B,C,n_delay)

xsize = length(A);
usize = size(B,2);
ysize = size(C,1);
dsize = 2; % number of disturbances

%% Define augmented system
% delay of n_delay time steps in 2nd component of u
% constant output disturbance

Adelay1 = [zeros(n_delay(1)-1,1), eye(n_delay(1)-1); zeros(1,n_delay(1))]; % delay component of A
Adelay2 = [zeros(n_delay(2)-1,1), eye(n_delay(2)-1); zeros(1,n_delay(2))]; % delay component of A
Adist = eye(dsize); % disturbance component of A

if (n_delay(1)==0)
    Aaug = [A, B(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(2),xsize), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(2)), Adist];

    Baug = [B(:,1), zeros(xsize,1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1);1];
        zeros(dsize,2)];
else
    Aaug = [A, B(:,1), zeros(xsize,n_delay(1)-1), B(:,2), zeros(xsize,n_delay(2)-1), zeros(xsize,dsize);
        zeros(n_delay(1),xsize), Adelay1, zeros(n_delay(1),n_delay(2)), zeros(n_delay(1),dsize);
        zeros(n_delay(2),xsize), zeros(n_delay(2),n_delay(1)), Adelay2, zeros(n_delay(2),dsize);
        zeros(dsize,xsize), zeros(dsize,n_delay(1)), zeros(dsize,n_delay(2)), Adist];

    Baug = [zeros(xsize,usize);
        [zeros(n_delay(1)-1,1); 1], zeros(n_delay(1),1);
        zeros(n_delay(2),1), [zeros(n_delay(2)-1,1); 1];
        zeros(dsize,usize)];
end

% disturbance to output matrix
Cd = eye(ysize,dsize); 

Caug = [C, zeros(ysize,sum(n_delay)), Cd];

Daug = 0;

end
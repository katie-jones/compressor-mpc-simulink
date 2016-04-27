[~,~,~,p,m] = const_mpc();

UW = [4e3 6e4];
YW = [1 1 0.1 5e3];

UWT = kron(eye(m),diag([UW,UW]'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
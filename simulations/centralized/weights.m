[~,~,~,p,m] = const_mpc();

UW = [2e4 2e5];
YW = [1 1 0.1 5e2];

UWT = kron(eye(m),diag([UW,UW]'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
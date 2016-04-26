[~,~,~,p,m] = const_mpc();

UW = [1.2e4 1.2e5];
YW = [1 0 8e3];

UWT = kron(eye(m),diag([UW]'));
YWT = kron(eye(p),diag(YW'));
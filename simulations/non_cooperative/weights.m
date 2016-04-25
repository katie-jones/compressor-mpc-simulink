[~,~,~,p,m] = const_mpc();

UW = [3e2 1e4];
YW = [1 0 1e6];

UWT = kron(eye(m),diag([UW]'));
YWT = kron(eye(p),diag(YW'));
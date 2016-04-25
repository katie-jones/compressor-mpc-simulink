[~,~,~,p,m] = const_mpc();

UW = [1e4 2e5];
YW = [1 1 1e-5 1e4];

UWT = kron(eye(m),diag(UW));
YWT = kron(eye(p),diag(YW'));

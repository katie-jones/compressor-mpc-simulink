[~,~,~,p,m] = const_mpc();

UW = [4e4 4e5];
YW = [1 1 1e-5 3e2];

UWT = kron(eye(m),diag(UW));
YWT = kron(eye(p),diag(YW'));

[~,~,~,p,m] = const_mpc();

UW = 5e2*[1e2 1e3];
YW = [200 1 1000 5];

UWT = kron(eye(m),diag(UW));
YWT = kron(eye(p),diag(YW'));

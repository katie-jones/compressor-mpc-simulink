[~,~,~,p,m] = const_mpc();

UW = 1e4*[1e2 1e3];
YW = [1000 1];

UWT = kron(eye(m),diag([UW]'));
YWT = kron(eye(p),diag(YW'));
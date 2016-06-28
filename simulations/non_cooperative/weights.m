[~,~,~,p,m] = const_mpc();

UW = 1e4*[1e2 1e3];
YW1 = [1e3 10];
YW2 = [1e4 1e2];

UWT = kron(eye(m),diag([UW]'));
YWT1 = kron(eye(p),diag(YW1'));
YWT2 = kron(eye(p),diag(YW2'));
[~,~,~,p,m] = const_mpc();

UW = 1e5*[1e2 1e3];
YW1 = [1000 10];
YW2 = [1000 10];

UWT = kron(eye(m),diag([UW]'));
YWT1 = kron(eye(p),diag(YW1'));
YWT2 = kron(eye(p),diag(YW2'));
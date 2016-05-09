UW = [1e2 1e4 1e2 1e4];
YW = [0 1 0 1];

UWT = kron(eye(m),diag(UW'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
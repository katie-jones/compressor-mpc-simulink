UW = 1e0*[1e1 1e6 1e1 1e6];
YW = [0 1e-4 1 1e-4];

UWT = kron(eye(m),diag(UW'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
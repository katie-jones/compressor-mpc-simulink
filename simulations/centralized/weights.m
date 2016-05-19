UW = 5e2*[1e1 6e2 1e1 6e2];
YW = [15 1 23 1];

UWT = kron(eye(m),diag(UW'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
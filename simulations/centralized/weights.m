UW = 1e2*[1e1 1e3 1e1 1e3];
YW = [5 1 10 1];

UWT = kron(eye(m),diag(UW'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
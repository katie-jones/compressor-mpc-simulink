UW = 2.5e2*[1e2 1e3 1e2 1e3];
YW = [200 1 1000 5];

UWT = kron(eye(m),diag(UW'));
% UWT = kron([1, 0; 0, 60],diag([UW,UW]));
YWT = kron(eye(p),diag(YW'));
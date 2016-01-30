% Various methods for computing an inverse of the Hessian

T = 100;
D = spdiags([-ones(T,1) ones(T,1)], [-1 0], T, T-1);
w = ceil(rand(T,1)*10);

idx = randperm(T-1,10);

% Fast way
a = 1./w(1:end-1) + 1./w(2:end);
b = 1./w(2:end-1);
cslogb = cumsum(log(b));
tic; V1 = tri_inverse2(a, -b, cslogb, idx); toc

% Reasonably fast way
A = D'*diag(1./w)*D;
R = sparse(chol(A));
tic;
V2 = R \ (R' \ E);
V2 = V2(idx,:);
toc;

% Slow way
V3 = inv(A);
V3 = V3(idx,idx);

% Also fast?
[rd,re] = tri_factor2(a,b);
E = eye(T-1);
tic; V4 = tri_solve2(rd,re,E(idx,:)); toc

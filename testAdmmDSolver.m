
randn('seed',1);
rand('seed',1);

m = 4;
n = 3;

% Differences on lattice
A = zeros(m,n); A(:,1:end-1) = 1;
B = zeros(m,n); B(:,2:end) = -1;
ai = find(A); bi = find(B);
D1 = [sparse(1:length(ai), ai, 1, length(ai), m*n) + ...
      sparse(1:length(bi), bi, -1, length(bi), m*n)];

A = zeros(m,n); A(1:end-1,:) = 1;
B = zeros(m,n); B(2:end,:) = -1;
ai = find(A); bi = find(B);
D2 = [sparse(1:length(ai), ai, 1, length(ai), m*n) + ...
      sparse(1:length(bi), bi, -1, length(bi), m*n)];

D = [D1; D2];
Y = randn(m*n,2);

opts = options('adm');
opts.maxit = 10;
lambda = 0.1;
out = admmDSolver(Y,D,lambda,opts);

cvx_begin
  variable X(m*n,2);
  minimize( 0.5*sum(sum_square(X-Y)) + lambda*sum(norms(D*X,[],2)) )
cvx_end

X0 = X;
X1 = out.X;
f0 = 0.5*sum(sum_square(X0-Y)) + lambda*sum(norms(D*X0,[],2));
f1 = 0.5*sum(sum_square(X1-Y)) + lambda*sum(norms(D*X1,[],2));

assert((f1-f0)/f0 < 1e-2);
disp('PASSED');
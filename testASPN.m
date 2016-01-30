
rand('seed',1);
T = 3000;
n = 10;
Y = rand(n,T)*2-1;

lambda = 10;
lam = lambda*ones(T-1,1);
opt = options('aspn');
out = runMethod(Y,lam,opt);

cvx_begin
  variable X(n,T)
  minimize( 0.5*sum(sum_square(X-Y)) + ...
            lambda*sum(norms(X(:,2:end)-X(:,1:end-1))) )
cvx_end

X0 = X;
X1 = out.X;
f0 = 0.5*sum(sum_square(X0-Y)) + lambda*sum(norms(X0(:,2:end)-X0(:,1:end-1)));
f1 = 0.5*sum(sum_square(X1-Y)) + lambda*sum(norms(X1(:,2:end)-X1(:,1:end-1)));
assert(abs(f1-f0)/f0 < 1e-8);

disp('PASSED');

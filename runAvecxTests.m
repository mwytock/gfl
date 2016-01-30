% Run 1/2||Avec(X)-y||_2^2 + lambda*TV(X) applications

randn('seed',1);
rand('seed',1);

problems = {problemTimeSeries(1000,20,[0.1 0.15],0.05)
            problemTimeSeries(10000,20,[0.1 0.1],0.05)};
lambdas = [1 10];
rho1 = [0.1 0.1];
rho2 = [10 10];

outs = cell(length(problems),1);
regopts = cell(length(problems),1);
opts = cell(length(problems),1);
cvx_X = cell(length(problems),1);
min_f = inf(length(problems),1);

for i=1:length(problems)
  y = problems{i}.y;
  A = problems{i}.A;
  T = length(y);
  n = size(A,2)/T;
  lam = lambdas(i)*ones(T-1,1);

  fprintf('T=%d\n');

  % ADMM w/ TV regularization
  regopts{i} = options('aspn');
  regopts{i}.printEvery = 1e10;
  regopts{i}.max_iters = 100;
  regopts{i}.k_max = 50;
  regopts{i}.stats = false;
  regopts{i}.verbose = 0;
  opts{i}{1} = options('adm');
  opts{i}{1}.rho = rho1(i);
  opts{i}{1}.reltol = 1e-6;
  opts{i}{1}.abstol = 1e-6;
  outs{i}{1} = runAvecx(A,y,lam,opts{i}{1},regopts{i});

  % vanilla ADMM
  e = ones(T*n,1);
  F = spdiags([-e e], [0 n], T*n, T*n);
  part = n*ones(T-1,1);
  w = ones(T,1);
  opts{i}{2} = options('adm');
  opts{i}{2}.rho = rho2(i);
  opts{i}{2}.max_iter = 1000;
  opts{i}{2}.reltol = 1e-6;
  opts{i}{2}.abstol = 1e-6;
  outs{i}{2} = admmAvecxSolver(A,F,y,w,lam(1),part,opts{i}{2});

  % CVX / SDPT3 low precision
  cvx_begin
    cvx_precision low
    variable X(n,T)
    minimize( 0.5*sum_square(y-A*vec(X)) +  ...
              lam(1)*sum(norms(X(:,2:end)-X(:,1:end-1))) )
  cvx_end
  outs{i}{3}.X = X;
  outs{i}{3}.obj = cvx_optval;

  % CVX / SDPT3 default precision
  cvx_begin
    cvx_precision default
    variable X(n,T)
    minimize( 0.5*sum_square(y-A*vec(X)) +  ...
              lam(1)*sum(norms(X(:,2:end)-X(:,1:end-1))) )
  cvx_end
  outs{i}{4}.X = X;
  outs{i}{4}.obj = cvx_optval;
  min_f(i) = cvx_optval;
end


save('avecx_tests', 'problems', 'lambdas', 'regopts', 'opts', 'outs');
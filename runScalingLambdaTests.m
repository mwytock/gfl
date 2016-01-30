% Run tests showing how ASPN scales as a function of problem size and
% sparsity

rand('seed', 1)

% Scaling as a function of sparsity (lambda)
T = 10000;
n = 10;
lambdas = logspace(log10(2),2,100);
Y = rand(n,T)*2-1;
k_max = 2500;

outs = cell(length(lambdas),1);
opts = cell(length(lambdas),1);
for i=1:length(lambdas)
  lambda = lambdas(i);
  fprintf('lambda=%f k_max=%d\n', lambda, k_max);
  lam = lambda*ones(T-1,1);
  opts{i} = options('aspn');
  opts{i}.verbose = 1;
  opts{i}.k_max = k_max;
  outs{i} = runMethod(Y,lam,opts{i});
  k_max = length(find(outs{i}.z>0));
end

nnzs = zeros(length(outs),1);
time = zeros(length(outs),1);
for i=1:length(outs)
  time(i) = outs{i}.time(end);
  nnzs(i) = length(find(outs{i}.z>0));
end

save('scaling_lambda', 'nnzs', 'time');

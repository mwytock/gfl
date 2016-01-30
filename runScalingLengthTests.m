% Run tests showing how ASPN scales as a function of problem size and
% sparsity

% Scaling as a function of sparsity (lambda)

Ts = round(logspace(2,6,100));
n = 10;
k = 10;

opt = options('aspn');
opt.k_max = 100;

outs = cell(length(Ts),1);
for i=1:length(Ts)
  T = Ts(i);
  Y = problemMultiChange(T,n,k);
  fprintf('T=%d\n', T);
  lam = sqrt(T)*ones(T-1,1);
  outs{i} = runMethod(Y,lam,opt);
end
save('scaling_length', 'Ts', 'n', 'k', 'opt', 'outs');

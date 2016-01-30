% Run 1/2||AX-Y||_F^2 + lambda*TV(X) applications

randn('seed',1);
rand('seed',1);

problems = {problemTimeSeries(1000,5)};
lambdas = [1];

algs = {'adm'};

outs = cell(length(problems),1);
regopts = cell(length(problems),1);
opts = cell(length(problems),1);
min_f = inf(length(problems),1);

for i=1:length(problems)
  Y = problems{i}.Y;
  A = problems{i}.A;
  [n T] = size(Y);
  p = size(A,2);
  lam = lambdas(i)*ones(T-1,1);
  for j=1:length(algs)
    fprintf('T=%d alg=%s\n', T, algs{j});
    regopts{i}{j} = options('pn');
    regopts{i}{j}.verbose = 0;
    opts{i}{j} = options(algs{j});
    outs{i}{j} = runAx(A,Y,lam,opts{i}{j},regopts{i}{j});
    %min_f(i) = min(min_f(i),outs{i}{j}.obj(end));
  end
end

save('ax_tests', 'lambdas', 'regopts', 'opts', 'outs');
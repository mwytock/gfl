% Run timing tests on real and synthetic data

randn('seed',1);
rand('seed',1);

problems = {problemMultiChange(1e3,1e3/100,1e3/1000)
            problemMultiChange(1e4,1e4/100,1e4/1000)
            problemBladder()
            problemLung()};
problems{5} = problems{3};
problems{6} = problems{4};

lambdas = [2 20 100 1000 2 150];
k_max = [50 200 20 20 700 300];

algs = {'aspn', 'pn', 'gflasso', 'pg', 'og', 'dr', 'lb', 'plb', 'sbb'};

outs = cell(length(problems),1);
opts = cell(length(problems),1);
min_f = inf(length(problems),1);
min_gap = inf(length(problems),1);

for i=1:length(problems)
  Y = problems{i};
  [n T] = size(Y);
  lam = lambdas(i)*ones(T-1,1);
  for j=1:length(algs)
    fprintf('T=%d alg=%s\n', T, algs{j});
    opts{i}{j} = options(algs{j});
    opts{i}{j}.k_max = k_max(i);
    outs{i}{j} = runMethod(Y,lam,opts{i}{j});
    min_f(i) = min(min_f(i),outs{i}{j}.obj(end));
    if isfield(outs{i}{j}, 'gap')
      min_gap(i) = min(min_gap(i),outs{i}{j}.gap(end));
    end
  end

  outi = outs{i};
  optsi = opts{i};
  save(['tests_' num2str(i)]);
end

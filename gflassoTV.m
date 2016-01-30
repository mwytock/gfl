function out = gflassoTV(Y, lam, opts)
opts.weights = ones(length(lam),1);
out = gflasso(Y',lam(1),opts);

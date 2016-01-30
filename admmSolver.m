function out=admmSolver(FX,REG,X0,opts)
t_start = tic;
Z = X0;
U = zeros(size(X0));

if opts.verbose >= 1
  fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
          'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

for i=1:opts.maxit
  X = FX.prox(Z-U);
  Zold = Z;
  Z = REG.prox(X+U);

  U = U+X-Z;
  out.obj(i) = FX.fgx(Z);
  out.time(i) = toc(t_start);

  r_norm = norm(X-Z, 'fro');
  s_norm = norm(opts.rho*(Z - Zold),'fro');
  eps_pri = sqrt(prod(size(X0)))*opts.abstol + ...
            opts.reltol*max(norm(X,'fro'), norm(Z,'fro'));
  eps_dual = sqrt(prod(size(X0)))*opts.abstol + ...
      opts.reltol*norm(opts.rho*U,'fro');

  if opts.verbose >= 1
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', i, ...
            r_norm, eps_pri, s_norm, eps_dual, out.obj(i));
  end

  if (r_norm < eps_pri && s_norm < eps_dual)
    break
  end
end

out.X = Z;
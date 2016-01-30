function out = admmDSolver(Y,D,lambda,opts)
[T n] = size(Y);

out.time0 = cputime;
R = chol(opts.rho*D'*D+speye(T));
Z = D*Y;
U = zeros(size(Z));
N = sqrt(prod(size(Z)));

for i=1:opts.maxit
  % X update
  X = R \ (R' \ (Y + opts.rho*D'*(Z-U))); 

  % Z update
  Zold = Z;

  Z = D*X+U;
  nz = sqrt(sum(Z.^2,2));
  Z = bsxfun(@times, pos(1-(lambda/opts.rho)./nz), Z);

  % U update
  U = U + (D*X-Z);

  r_norm = norm(D*X-Z,'fro');
  s_norm = norm(opts.rho*(Z-Zold),'fro');
  eps_pri = opts.abstol*N + opts.reltol*max(norm(D*X,'fro'), norm(Z,'fro'));
  eps_dual = opts.abstol*N + opts.reltol*norm(opts.rho*U,'fro');
  out.f(i) = 0.5*norm(X-Y,'fro')^2 + lambda*sum(sqrt(sum((D*X).^2,2)));
  out.time(i) = cputime;

  if opts.verbose >= 1
      fprintf('%d\t%f\t%f\t%f\t%f\t%f\n', i, r_norm, eps_pri, ...
              s_norm, eps_dual, out.f(i));
  end

  if r_norm < eps_pri && s_norm < eps_dual
    break
  end
end

out.X = X;
out.Z = Z;

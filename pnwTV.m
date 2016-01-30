function out = pnwTV(Y, w, lam, opt)

t_start = tic;
[n T] = size(Y);
D = spdiags([-ones(T,1) ones(T,1)], [-1 0], T, T-1);
YD = Y*D;

a = 1./w(1:end-1) + 1./w(2:end);
b = 1./w(2:end-1);
cslogb = cumsum(log(b));
E = eye(T-1);

% Initialize
if (isfield(opt, 'z0'))
  z = opt.z0;
else
  z = zeros(T-1,1);
end

[rd,re] = tri_factor2(a+z,b);
U = tri_solve2(rd,re,YD);
f = 0.5*sum(sum(YD.*U)) + z'*lam.^2/2;
num_ls = 0;

for i=1:opt.max_iters
  g = -0.5*sum(U.^2)' + lam.^2/2;
  idx = find(~(z == 0 & g > opt.eps));
  out.gzI(i) = norm(g(idx));

  if opt.verbose >= 1
    fprintf('%4d\t%10f\t%6e\t%10d\n', i, f, out.gzI(i), length(idx));
  end

  if out.gzI(i) < opt.tol
    break
  end

  V = tri_solve2(rd,re,E(idx,:));
  H = V(:,idx) .* (U(:,idx)'*U(:,idx));
  d = -H\g(idx);

  alpha = 1;
  z_next = z;
  df = g(idx)'*d;
  for j=1:opt.max_ls_iters
    num_ls = num_ls + 1;
    z_next(idx) = z(idx) + alpha*d;
    z_next(z_next < 0) = 0;
    [rd,re] = tri_factor2(a+z_next,b);
    U_next = tri_solve2(rd,re,YD);
    f_next = 0.5*sum(sum(YD.*U_next)) + z_next'*lam.^2/2;

    if opt.verbose >= 2
      fprintf('  LS %4d\t%10f\t%10f\n', j, alpha, f_next);
    end

    if f_next <= f + 0.1*alpha*df
      break
    end

    if opt.ls_quadratic && j == 1
      alpha = -df*alpha^2/(2*(f_next - f - df*alpha));
    else
      alpha = alpha*opt.beta;
    end
  end

  % TODO(mwytock): Fix the stopping criterion to use duality gap or similar
  if (f - f_next < opt.tol)
    break
  end

  z = z_next;
  U = U_next;
  f = f_next;
end

if i==opt.max_iters
  save('pnwTV_fail', 'Y', 'w', 'lam', 'opt');
  error('pnwTV hit iteration limit');
end

out.z = z;
out.f = f;
out.X = Y - (U*D').*repmat(1./w',n,1);

if opt.verbose >= 1
  fprintf('Total number of line searches: %d\n', num_ls);
end

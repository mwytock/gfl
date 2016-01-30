function out = pnTV(Y, lam, opt)

t_start = cputime;
[n T] = size(Y);
D = spdiags([-ones(T,1) ones(T,1)], [-1 0], T, T-1);

if opt.verbose >= 1
  fprintf('%4s\t%10s\t%10s\t%10s\n', 'iter', 'obj', 'gap', '#free');
end

% Initialize
if (isfield(opt, 'z0'))
  z = opt.z0;
else
  z = zeros(T-1,1);
end

if length(lam) ~= size(Y,2)-1
   lam=lam*ones(size(Y,2)-1,1);
end

num_ls = 0;
[rd,re] = dtd_factor(z);
U = dtd_solve(rd,re,Y*D);
f = 0.5*sum(sum(Y.*(U*D'))) + z'*lam.^2/2;

for i=1:opt.max_iters
  g = -0.5*sum(U.^2)' + lam.^2/2;

  idx0 = find(z > 0);
  idx1 = find(z == 0 & g < opt.eps);
  idx = [idx0; idx1];
  if isfield(opt, 'k_max') && length(idx) > opt.k_max
    [val, ind] = sort(abs(g(idx1)), 'descend');
    idx = [idx0; idx1(ind(1:opt.k_max-length(idx0)))];
  end

  out.dd_obj(i) = f;
  out.iterTime(i) = cputime;
  [out.obj(i), out.gap(i)] = primal_objval2(Y,U,lam(1));
  out.statsTime(i) = cputime;

  if cputime-t_start > opt.max_time
    break
  end

  if opt.verbose >= 1
    fprintf('%4d\t%10f\t%6e\t%10d\n', i, f, out.gap(i), length(idx));
  end

  if out.gap(i) < opt.tol
    break
  end

  idx = sort(idx);
  V = dtd_inverse(z, idx);
  H = V .* (U(:,idx)'*U(:,idx));
  H0 = spdiags(spdiags(inv(H),[-10:10]), [-10:10], size(H,1), size(H,2));
  d = -H\g(idx);
  %d = -H0*g(idx);
  %disp(num2str(norm(inv(H0)-inv(H),'fro')))
  

  alpha = 1;
  z_next = z;
  df = g(idx)'*d;
  for j=1:opt.max_ls_iters
    num_ls = num_ls + 1;
    z_next(idx) = z(idx) + alpha*d;
    z_next(z_next < 0) = 0;
    [rd,re] = dtd_factor(z_next);
    U_next = dtd_solve(rd,re,Y*D);
    f_next = 0.5*sum(sum(Y.*(U_next*D'))) + z_next'*lam.^2/2;

    if opt.verbose >= 2
      fprintf('  LS %4d\t%10f\t%10f\n', j, alpha, f_next);
    end

    if f_next <= f + 0.1*alpha*df
      break
    end

    if opt.ls_quadratic && j==1
      alpha = -df*alpha^2/(2*(f_next - f - df*alpha));
    else
      alpha = alpha*opt.beta;
    end
  end

  z = z_next;
  U = U_next;
  f = f_next;
end

V = projInfTwoBall(lam(1),U);
out.X = Y-rightMultByEtrans(V);
out.z = z;
out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;

if opt.verbose >= 1
  fprintf('Total number of line searches: %d\n', num_ls);
end

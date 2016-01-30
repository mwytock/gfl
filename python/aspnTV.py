function out = aspnTV(Y, lam, opt)

t_start = cputime;
[n T] = size(Y);
YD = Y(:,1:end-1)-Y(:,2:end);
Ycs = [zeros(n,1) cumsum(Y,2)];
Y2cs = [zeros(n,1) cumsum(Y.^2,2)];

% Initialize
if (isfield(opt, 'z0'))
  z = opt.z0;
else
  z = zeros(T-1,1);
end

if opt.verbose >= 1
  fprintf('%4s\t%10s\t%10s\t%10s\n', 'iter', 'obj', 'gz', '#free');
end

f = sum(Y2cs(:,end)) - sum(Ycs(:,end).^2/T);
out2.f = 0;
z = zeros(T-1,1);

for i=1:opt.max_iters,
  % compute entries of working set
  [rd,re] = dtd_factor(z);
  U = dtd_solve(rd,re,YD);
  g = -0.5*sum(U.^2)' + lam.^2/2;

  nZ = length(find(z>0));
  idx = find(~(z == 0 & g > opt.eps));
  if isfield(opt, 'k_max') && length(idx) > opt.k_max
    [val, ind] = sort(abs(g(idx)), 'descend');
    idx = unique([idx(ind(1:opt.k_max-nZ)); find(z>0)]);
  end

  % compute function
  out.dd_obj(i) = out2.f + 0.5*sum(f);
  out.gz(i) = norm(g(idx));
  out.iterTime(i) = cputime;
  [out.obj(i), out.gap(i)] = primal_objval2(Y,U,lam(1));
  out.statsTime(i) = cputime;

  if opt.verbose >= 1
    fprintf('%4d\t%6e\t%6e\t%5d %5d\n', i, out.obj(i), out.gap(i), nZ, length(idx));
  end

  if cputime-t_start > opt.max_time
    break
  end

  if out.gz(i) < opt.tol
    break;
  end

  if i > 1 && abs(out.obj(i) - out.obj(i-1)) < opt.tol
    break;
  end

  % compute reduced problem
  t = [0 idx' T];
  T0 = length(idx)+1;
  Y0 = zeros(n,T0);
  w = zeros(T0,1);
  f = zeros(T0,1);
  for j=1:T0,
    w(j) = t(j+1) - t(j);
    Y0(:,j) = (Ycs(:,t(j+1)+1) - Ycs(:,t(j)+1))/w(j);
    f(j) = sum(Y2cs(:,t(j+1)+1) - Y2cs(:,t(j)+1)) - sum(Y0(:,j).^2*w(j));
  end

  % solve reduced problem
  opt2 = opt;
  opt2.verbose = 2;
  opt2.z0 = z(idx);
  opt2.tol = opt.tol/10;
  out2 = pnwTV(Y0,w,lam(idx),opt2);
  z(idx) = out2.z;
end
out.z = z;
out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;

V = projInfTwoBall(lam(1),U);
X = Y - rightMultByEtrans(V);
out.X = X;
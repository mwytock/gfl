function history = admmAvecxSolver(A, F, b, w, lambda, part, params)
t_start = tic;
rho = params.rho;

%% Global constants and defaults

if ~isfield(params, 'quiet') params.quiet = false; end
if ~isfield(params, 'max_iter') params.max_iter = 1000; end
if ~isfield(params, 'abstol') params.abstol = 1e-4; end
if ~isfield(params, 'reltol') params.reltol = 1e-2; end

abstol = params.abstol;
reltol = params.reltol;

%% Data preprocessing

[m,n] = size(A);

% save a matrix-vector multiply
Atb = A'*spdiags(w,0,m,m)*b;
if sum(part) > n
  error('invalid partition');
end

% cumulative partition
cum_part = cumsum(part);

%% ADMM solver

x = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

% pre-factor
L = chol(A'*spdiags(w,0,m,m)*A + rho*F'*F, 'lower');
U = L';

if ~params.quiet
  fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
          'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

for k=1:params.max_iter
  % x-update
  q = Atb + rho*F'*(z - u);    % temporary value
  x = U \ (L \ q);

  % z-update
  zold = z;
  start_ind = 1;
  Fx = F*x;
  z = Fx;

  % Hacks to vectorize
  T = m;
  UFX = reshape(u+Fx,[],T);
  FX = reshape(Fx,[],T);
  normsUFX = norms(UFX);
  kappa = lambda/rho;
  Z = zeros(n/m, m);
  Z(:,1:T-1) = pos(1 - kappa./repmat(normsUFX(:,1:T-1),n/m,1)).*UFX(:,1:T-1);
  Z(:,T) = FX(:,T);
  z = reshape(Z,[],1);

  u = u + (Fx - z);

  % diagnostics, reporting, termination checks
  history.time(k) = toc(t_start);
  history.objval(k)  = objective(A, F, b, w, lambda, cum_part, x, z);
  history.r_norm(k)  = norm(Fx - z);
  history.s_norm(k)  = norm(-rho*(z - zold));
  history.eps_pri(k) = (sqrt(n)*abstol + reltol*max(norm(x), norm(-z)));
  history.eps_dual(k)= (sqrt(n)*abstol + reltol*norm(rho*u));
  t_start = tic;

  if ~params.quiet
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
  end

  if (history.r_norm(k) < history.eps_pri(k) && ...
      history.s_norm(k) < history.eps_dual(k))
    break;
  end
end

history.time = cumsum(history.time);
history.x = x;

end

function z = shrinkage(x, kappa)
  z = pos(1 - kappa/norm(x))*x;
end

function p = objective(A, F, b, w, lambda, cum_part, x, z)
  obj = 0;
  start_ind = 1;
  z = F*x;
  for i = 1:length(cum_part),
    ii = start_ind:cum_part(i);
    obj = obj + norm(z(ii));
    start_ind = cum_part(i) + 1;
  end
  r = (A*x - b).*sqrt(w);
  p = ( 0.5*r'*r + lambda*obj );
end

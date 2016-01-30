function out=runAvecx(A,y,lam,opts,regopts)
% RUNFXRX -- solves min 1/2||A*vec(X)-y||_2^2 + TV(X)
  [T Tn] = size(A);
  n = Tn/T;

  % Initialize state
  global state;
  state.z0 = zeros(T-1,1);

  % AA^T is diagonal so we just invert it and get a block diagonal matrix
  % which we can use in O(n) time on each iteration
  scale = 1/opts.rho;
  state.B = speye(Tn) - A'*inv(speye(T)/scale + A*A')*A;
  state.Aty = A'*y*scale;

  REG.fval = @(X) lam(1)*onetwonorm(X(:,1:T-1)-X(:,2:T));
  REG.prox= @(Z) proxRx(Z,scale*lam,regopts);
  FX.fgx = @(X) objFunGrad(X, REG, A, y);
  FX.prox = @(Z) proxFx(Z);

  X0 = zeros(n,T);
  out=admmSolver(FX,REG,X0,opts);
end

function [f,g]=objFunGrad(X,REG,A,y)
  ax=A*vec(X)-y;
  f=0.5*norm(ax,'fro')^2 + REG.fval(X);
  if (nargout > 1)
    g = 0;  % Not used currently
  end
end

function X=proxFx(Z)
  global state;
  x = state.B*(state.Aty + vec(Z));
  X = reshape(x,size(Z));
end

function X=proxRx(Z,lam,regopts)
  global state;
  regopts.z0 = state.z0;
  out=runMethod(Z,lam,regopts);
  X=out.X;
  state.z0 = out.z;
end

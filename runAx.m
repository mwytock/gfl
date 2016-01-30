function out=runAx(A,Y,lam,opts,regopts)
% RUNFXRX -- solves min 1/2||AX-Y||_F^2 + TV(X)
  [n,p] = size(A);
  T = size(Y,2);

  REG.fval = @(X) lam(1)*onetwonorm(X(:,1:T-1)-X(:,2:T));
  REG.prox= @(Z,scale) proxFn(Z,scale*lam,regopts);

  % TODO(mwytock): Pre factor
  ATA = A'*A;
  ATY = A'*Y;
  FX.fgx = @(X) objFunGrad(X, REG, A, Y);
  FX.prox = @(Z,scale) (scale*ATA + eye(p))\(scale*ATY+Z);

  X0 = zeros(p,T);
  switch opts.alg
    case 'adm'
      out=admmSolver(FX,REG,X0,opts);
    case 'og'
      out=adaptiveFista(FX,REG,X0,opts);
    case 'trip'
      out=TRIP(FX,REG,X0,opts);
    otherwise
      error('runAx: Unrecognized solver requested');
  end
end

function [f,g]=objFunGrad(X,REG,A,Y)
  ax=A*X-Y;
  f=0.5*norm(ax,'fro')^2 + REG.fval(X);
  if (nargout > 1)
    g = A'*ax;
  end
end

function X=proxFn(Z,lam,regopts)
  out=runMethod(Z,lam,regopts);
  X=out.X;
end

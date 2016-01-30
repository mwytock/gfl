function out=runFxRx(alg, x0, lam, varargin)
% RUNFXRX -- solves min f(x) + reg(x)
   
   opts=options('pn');
   opts.verbose=0;
   n = size(x0,2);
   REG.fval = @(X) lam*onetwonorm(X(:,1:n-1)-X(:,2:n));
   REG.prox= @(X,al) proxFn(X,al*lam,opts);
   FX.fgx = @(X) objFunGrad(X, REG, varargin{:});
   FX.prox = @(X) objFunProx(X, varargin{:});

   switch alg
     case 'adm'
       out=admmSolver(FX, REG,solopt);
     case 'og'
       out=adaptiveFista(FX,REG,solopt);
     case 'trip'
       solopt=options('trip');
       out=TRIP(FX,REG,x0,solopt);
     otherwise 
       error('runFxRx: Unrecognized solver requested');
   end
end

function [f,g]=objFunGrad(X, REG, varargin)
   A=varargin{1};
   Y=varargin{2};
   ax=A*X-Y;
   f=0.5*norm(ax,'fro')^2 + REG.fval(X);
   if (nargout > 1)
      g = A'*ax;
   end
end

function X=proxFn(X,lam,opts)
   out=pnTV(X,lam,opts);
   X=out.X;
end

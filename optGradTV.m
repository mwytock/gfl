function out=optGradTV(Y,lam,opts)
  t_start = cputime;
% OPTGRADTV -- run accelerated gradient for GFL-TV op
   UET = @(X) rightMultByEtrans(X);
   U = initU(Y,lam,opts);               % to allow unified initialization
   [d,n]=size(Y);
   % TODO(mwytock): Fix to take lambda as vector
   lam = lam(1);
   Z=U;

   for i=1:opts.maxit
      th=2/(i+1);
      Up=U;
      V=(1-th)*U+th*Z;
      X = UET(V)-Y;
      XE=X(:,2:n)-X(:,1:n-1);
      U = projInfTwoBall(lam,V+0.25*XE);
      Z = Up + (1/th)*(U-Up);
      out.iterTime(i) = cputime;
      [out.obj(i) out.gap(i)] = primal_objval2(Y,U,lam);
      out.statsTime(i) = cputime;
      if cputime-t_start > opts.max_time
        break
      end
   end
   out.U=U;
   out.X=-X;
   out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;
   out.gap(1) = inf;
end

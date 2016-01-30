function out=drTV(Y,lam,opts)
% DRTV -- run projected gradient for GFL-TV op
   t_start = cputime;
   UET = @(X) rightMultByEtrans(X);
   U = initU(Y,lam,opts);               % to allow unified initialization
   [d,n]=size(Y);
   YE= Y(:,1:n-1)-Y(:,2:n);             % diff operator on cols of Y

   % TODO(mwytock): Fix to take lambda as vector
   lam = lam(1);

   [rd,re]=tri_factor(ones(n-1,1));

   for i=1:opts.maxit
      R1=2*projInfTwoBall(lam,U)-U;
      R2=R1+YE;
      tri_solve(rd,re,R2);
      R2 = 2*R2-R1;
      U  = 0.5*(R2+U);
      out.iterTime(i) = cputime;
      [out.obj(i) out.gap(i)] = primal_objval2(Y,U,lam);
      out.statsTime(i) = cputime;
      if cputime-t_start > opts.max_time
        break
      end
   end
   out.U=projInfTwoBall(lam,U);
   out.U=U;
   out.X=Y-UET(U);
   out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;
end
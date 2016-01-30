function U=initU(Y,lam,opts)

   [d,n]=size(Y);
   if strcmp(opts.alg, 'spg')
      U=Y(:,1:n-1);
   else
      YE= Y(:,1:n-1)-Y(:,2:n);
      [rd,re]=tri_factor(zeros(n-1,1));
      U=YE;
      tri_solve(rd,re,U);
   end
end

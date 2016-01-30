function out=adaptiveFista(FX,REG,opts)

   Z=X;                                 % TODO INIT
   for i=1:opts.maxit
      th=2/(i+1);
      Xp=X;
      U=(1-th)*X+th*Z;
      al=getStepsize(U,X,th);
      [f,g]=FX.fgx(X);
      X = REG.prox(U - al*G);
      Z = Up + (1/th)*(U-Up);
      out.obj(i)=f;
   end
   out.X=X;
end
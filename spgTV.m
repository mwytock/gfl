function out=spgTV(Y,lam,opts)
% SPGTV -- run spectral projected gradient for GFL-TV op

   UET = @(X) rightMultByEtrans(X);

   U = initU(Y,lam,opts);               % to allow unified initialization

   % TODO(mwytock): Fix to take lambda as vector
   lam = lam(1);
   [d,n]=size(Y);

   prx  = @(X) projspg(X,lam,d,n-1);
   gfx  = @(X) spgGFX(X,d,n-1,Y);
   if opts.stats
      out.obj(1)=0.5*norm(UET(U)-Y,'fro')^2;
   end
   % hack cleanSPG to maintain stats only if asked to...
   [U,f,out.pr,out.it,~,info]=cleanSPG(gfx,U(:),prx,opts.maxit);
   out.obj = [out.obj; info(:)];
   U=reshape(U,d,n-1);
   out.U=U;
end

function [f,g]=spgGFX(U,d,m,Y)
   U=reshape(U,d,m);
   T=rightMultByEtrans(U)-Y;
   f=0.5*norm(T,'fro')^2;
   if (nargout > 1)
      n=size(T,2);
      g = T(:,1:n-1)-T(:,2:n);
      g = g(:);
   end
end

function U=projspg(V,lam,d,m)
   V=reshape(V,d,m);
   nv=sum(V.^2,1);nv=sqrt(nv);
   nv=max(lam,nv);
   U=bsxfun(@rdivide,V,nv/lam);
   U=U(:);
end

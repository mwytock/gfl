function out=zictv(Y,lam,mit,alg)

   addpath('./tri_solver');
   [d,n]=size(Y);
   e=ones(n,1);
   UET = @(x) uet(x);
   YE= Y(:,1:n-1)-Y(:,2:n);

   %out.E=E;
   %R=chol(EtE);
   %U = R \ (R' \ (YE'));
   %U=U';

   [rd,re]=tri_factor(zeros(n-1,1));
   U=YE;
   U(1,1)=YE(1,1);                      % CRUCIAL: TO FORCE MATLAB TO COPY
                                        % YE to U, else matlab
                                        % act. overwrites YE!!
   tri_solve(rd,re,U);


   itn=inftwonorm(U);
   if (itn <= lam)
      fprintf('|YE*inv(EtE)| = %E <= %E\n', itn,lam);
      out.U=U;
      return;
   else
      fprintf('|YE*inv(EtE)| = %E > %E\n', itn,lam);
   end

   switch alg
      case 'pg'
        for i=1:mit
           T = UET(U)-Y;                % X = -T, i.e., X = Y-UE'
           XE=T(:,2:n)-T(:,1:n-1);      % X*E
           U = proj(lam,U+0.25*XE);
           out.obj(i)=0.5*norm(T,'fro')^2;
           out.gap(i)=lam*onetwonorm(XE)-sum(sum(U.*XE));
        end
        f=0.5*norm(T,'fro')^2;
        out.U=U;
     case 'spg'
       prx  = @(x) projspg(x,lam,d,n-1);
       gfx  = @(x) spgGFX(x,d,n-1,Y);
       out.obj(1)=0.5*norm(UET(U)-Y,'fro')^2;
       U=Y(:,1:n-1);
       [U,f,out.pr,out.it,jnk,info]=cleanSPG(gfx,U(:),prx,mit);
       out.obj = [out.obj; info(:)];
       U=reshape(U,d,n-1);
       out.U=U;
     case 'og'
       Z=U;
       for i=1:mit
          th=2/(i+1);
          Up=U;
          V=(1-th)*U+th*Z;
          X = UET(V)-Y;
          XE=X(:,2:n)-X(:,1:n-1);
          U = proj(lam,V+0.25*XE);
          Z = Up + (1/th)*(U-Up);
          out.obj(i)=0.5*norm(X,'fro')^2;
          out.gap(i)=lam*onetwonorm(XE)-sum(sum(U.*XE));
       end
       f=0.5*norm(X,'fro')^2;
       out.U=U;
     case 'dr'
       %IE=chol(spdiags([-e 3*e -e], -1:1, n-1, n-1));
       %R1 = @(Z) refl1old(Z,YE,IE);

       [rd,re]=tri_factor(ones(n-1,1));

       %R1 = @(Z) refl1(Z,YE,rd,re);
       %R2 = @(Z) refl2(Z,lam);

       % made the calls inline
       for i=1:mit
          R1=2*proj(lam,U)-U;
          R2=R1+YE;
          tri_solve(rd,re,R2);
          R2 = 2*R2-R1;
          U  = 0.5*(R2+U);
          %U = 0.5*(R1(R2(U)) + U);
          if 0
             V=proj(lam,U);                % needed only statistics ---avoid
             T=UET(V)-Y;
             TE=T(:,1:n-1)-T(:,2:n);
             out.obj(i)=0.5*norm(T,'fro')^2;
             out.gap(i)=lam*onetwonorm(TE)+sum(sum((V.*TE)));
          end
       end
       out.U=proj(lam,U);
       f=0.5*norm(UET(U)-Y,'fro')^2;
     case 'lb'
       % run lbfgs-b on the dual of the dual
       N=(n-1);
       l=zeros(n-1,1);
       u=inf*ones(N,1);
       gfx = @(z) lbfgfx(z,Y,YE,lam);
       opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10, 'maxIts', mit, ...
                         'max_time', 1e4);
       [z, f, out.info] = lbfgsb(gfx, @(x) 1, l, u, opts);
       [rd,re]=tri_factor(z);
       U = YE;
       tri_solve(rd,re,U);
       out.U=U;
       out.z=z;
       out.obj=-out.info.err(:,1);
     case 'sr1'
       T = UET(U)-Y;
       prox = @(x) proj(lam,x);
       opts    = struct('nmax', mit);
       [xk,nit, errStruct, defaultOpts, stepsizes]=zeroSR1(fcn,grad,h,prox,opts);
     case 'ad'
       fprintf('ADMM: TODO\n');
       f=inf;
     case 'plb'
       gfx = @(z) lbfgfx(z,Y,YE,lam);
       U=Y(:,1:n-1);
       l=zeros(n-1,1);
       u=inf*ones(n-1,1);
       z0=zeros(n-1,1);
       opt=options();
       opt.maxit=mit;
       opt.verbose=0;
       out=plb(gfx,z0,l,u,opt);
       z=out.x;
       f=out.obj(end);
       [rd,re]=tri_factor(z);
       U = YE;
       tri_solve(rd,re,U);
       out.U=U;
       out.z=z;
     case 'bb'
       gfx = @(z) lbfgfx(z,Y,YE,lam);
       l=zeros(n-1,1);
       z0=zeros(n-1,1);
       opt=options();
       opt.maxit=mit;
       opt.verbose=1;
       out=runbb(gfx,z0,opt);
       z=out.x;
       f=out.obj(end);
       [rd,re]=tri_factor(z);
       U = YE;
       tri_solve(rd,re,U);
       out.U=U;
       out.z=z;
     case 'pn'
       % mainly just convert between suvrit's notation and ours
       params.verbose = 1;
       params.max_iters = mit;
       params.max_ls_iters = 10;
       params.eps = 1e-8;
       params.tol = 1e-4;
       params.beta = 0.5;
       params.ls_quadratic = 1;
       [z, f, U, history,X] = prox_mtv(Y', lam, zeros(size(Y,2)-1,1), params);
       U = -U';
       out.gap = history.gap;
       out.U = U;
       out.z = z;
       out.obj = history.objval;
       f=0.5*norm(UET(U)-Y,'fro')^2;
     otherwise
       fprintf('WARNING: Unknown method requested\n');
       f=inf;
   end
   out.X=Y-UET(U);
   fprintf('Obj=%d, |U| = %d\n',f,inftwonorm(out.U));
end

% maybe improve this to avoid so much memory traffic
function V=uet(U)
   n=size(U,2);
   K=U(:,2:n)-U(:,1:n-1);
   V=[U(:,1),K,-U(:,n)];
end

function [f,g]=lbfgfxold(z,Y,YE,EtE,lam)
   n=size(EtE,1);
   K=chol(spdiags(z, 0, n,n) + EtE);
   U = K \ (K' \ (YE'));
   U = U';

   nrms=sum(U.^2,1);nrms=nrms(:);
   l=lam*lam*ones(length(z),1);
   f=0.5*norm(uet(U)-Y,'fro')^2 + 0.5*sum(z.*(nrms-l));
   f=-f;
   if (nargout > 1)
      T=-0.5*U.*U; % Ki*YE'*YE*Ki;
      st=sum(T,1);
      g = 0.5*l + st';
      %g = -g;
   end
end


% SUVRIT: This version has a bug---the one without the fast tri_solve gives
% the right answer --- the one with seems to be wrong (definitely some kind
% of transposition error)
function [f,g]=lbfgfx(z,Y,YE,lam)
   [rd,re]=tri_factor(z);
   U=YE;
   U(1,1)=YE(1,1);                    % FORCE COPY
   tri_solve(rd,re,U);

   nrms=sum(U.^2,1);nrms=nrms(:);
   l=lam*lam*ones(length(z),1);
   f=0.5*norm(uet(U)-Y,'fro')^2 + 0.5*sum(z.*(nrms-l));
   f=-f;
   if (nargout > 1)
      T=-0.5*U.*U; % Ki*YE'*YE*Ki;
      st=sum(T,1);
      g = 0.5*l + st';
      %g = -g;
   end
end

function X=refl1old(Z,YE,IE)
   PX= IE \ (IE' \ (Z+YE)');
   PX=PX';
   X=2*PX-Z;
end

function X=refl1(Z,YE,rd,re)
% 2P-I, where P solves h|X-Z|^2+ h|XE'-Y|^2
   X=Z+YE;
   X(1,1)=Z(1,1);                       % FORCE COPY!!!!
   tri_solve(rd,re,X);
   X = 2*X-Z;
end

function X=refl2(Z,lam)
% 2P-I, where P is proj used below
   X=2*proj(lam,Z)-Z;
end

function [f,g]=spgGFX(U,d,m,Y)
   U=reshape(U,d,m);
   T=uet(U)-Y;
   f=0.5*norm(T,'fro')^2;
   if (nargout > 1)
      n=size(T,2);
      g = T(:,1:n-1)-T(:,2:n);
      g = g(:);
   end
end

function U=proj(lam,V)
   nv=sum(V.^2,1);nv=sqrt(nv);
   nv=max(lam,nv);
   U=bsxfun(@rdivide,V,nv/lam);
end

function U=projspg(V,lam,d,m)
   V=reshape(V,d,m);
   nv=sum(V.^2,1);nv=sqrt(nv);
   nv=max(lam,nv);
   U=bsxfun(@rdivide,V,nv/lam);
   U=U(:);
end

function d=inftwonorm(X)
   d=sum(X.^2,1);d=full(d);d=sqrt(d);
   d=max(d);
end

function d=onetwonorm(X)
   d=sum(X.^2,1);d=full(d);d=sqrt(d);
   d=sum(d);
end
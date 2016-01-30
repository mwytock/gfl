function out=plbTV(Y,lam,opts)
% PLBTV -- run plbfgs for GFL-TV op

   U = initU(Y,lam,opts);               % to allow unified initialization
   [d,n]=size(Y);
   YE= Y(:,1:n-1)-Y(:,2:n);             % diff operator on cols of Y

   % TODO(mwytock): Fix to take lambda as vector
   lam = lam(1);

   % run pbl on the dual of the dual
   gfx = @(z) lbfgfx(z,Y,YE,lam);
   fprimal = @(z) primal_objval(z,Y,YE,lam);
   l=zeros(n-1,1);
   u=inf*ones(n-1,1);
   z0=zeros(n-1,1);
   out=plb(gfx,fprimal,z0,l,u,opts);
   z=out.x;
   f=out.obj(end);
   [rd,re]=tri_factor(z);
   U = YE;
   tri_solve(rd,re,U);
   out.U=U;
   out.X=Y-rightMultByEtrans(U);
   out.z=z;
   out.dd_obj = out.obj;

   % do time accounting, remove statsTime
   out.obj = out.stats;
   out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - out.startTime;

   %fprintf('Primal: %E\t Dual=%E, |U| = %d\n',objval(out.X,Y,lam),f,inftwonorm(out.U));
end

function [f,g,info]=lbfgfx(z,Y,YE,lam)
   [rd,re]=tri_factor(z);
   U=YE;
   U(1,1)=YE(1,1);                    % FORCE COPY
   tri_solve(rd,re,U);

   nrms=sum(U.^2,1);nrms=nrms(:);
   l=lam*lam*ones(length(z),1);
   f=0.5*norm(rightMultByEtrans(U)-Y,'fro')^2 + 0.5*sum(z.*(nrms-l));
   f=-f;
   if (nargout > 1)
      T=-0.5*U.*U; % Ki*YE'*YE*Ki;
      st=sum(T,1);
      g = 0.5*l + st';
      %g = -g;
   end
 end

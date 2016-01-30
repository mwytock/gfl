function out=sbbTV(Y,lam,opts)
% SBBTV -- run SBB for GFL-TV op

   U = initU(Y,lam,opts);               % to allow unified initialization
   [d,n]=size(Y);
   YE= Y(:,1:n-1)-Y(:,2:n);             % diff operator on cols of Y

   % TODO(mwytock): Fix to take lambda as vector
   lam = lam(1);
   t_start = cputime;

   % run sbb on the dual of the dual
   N=n-1;
   l=zeros(N,1);
   u=inf*ones(N,1);
   gfx = @(z) lbfgfx(z,Y,YE,lam);
   fprimal = @(z) primal_objval(z,Y,YE,lam);
   lb=zictv(Y,lam,2,'lb');
   out=runbb(gfx,fprimal,lb.z,opts);
   z=out.x;
   f=out.obj(end);
   [rd,re]=tri_factor(z);
   U = YE;
   tri_solve(rd,re,U);
   out.U=U;
   out.X=Y-rightMultByEtrans(U);
   out.z=z;
   out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;
   %fprintf('Primal: %E\t Dual=%E, |U| = %d\n',objval(out.X,Y,lam),f,inftwonorm(out.U));
end

function [f,g]=lbfgfx(z,Y,YE,lam)
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

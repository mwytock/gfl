function [U t]=projInfTwoBall(lam,V)
% PROJINFTWOBALL -- project V onto inf-2-mixed norm ball
% 
% U=projInfTwoBall(lam,V)
%
% lam -- scalar for now
%
   tic;
   nv=sum(V.^2,1);nv=sqrt(nv);
   nv=max(lam,nv);
   U=bsxfun(@rdivide,V,nv/lam);
   t=toc;
end

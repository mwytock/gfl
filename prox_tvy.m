function X=prox_tvy(Y,prox)
% PROX_TVZ -- implements prox on z-slices
%
% Author: Suvrit Sra
%
% See also: prox_tvx, prox_mtv
%

   [nx,ny,nz]=size(Y);

   for i=1:ny                           % can be done in parallel
       Z=reshape(Y(:,i,:),nx,nz);
      X(:,i,:)=prox(Z);
   end
end
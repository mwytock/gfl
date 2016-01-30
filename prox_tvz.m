function X=prox_tvz(Y,prox)
% PROX_TVZ -- implements prox on z-slices
%
% Author: Suvrit Sra
%
% See also: prox_tvy, prox_mtv
%

   [nx,ny,nz]=size(Y);

   for i=1:nz                           % can be done in parallel
      X(:,:,i)=prox(Y(:,:,i));
   end
end
function out=denoise3d(Y,lam,mu,opts)
% DENOISE3D --- denoise a 3D Tensor
%
% See also: prox_tvz, prox_tvy
%
%
   mit=20;
   proxz = @(X) proxZ(X,lam,opts);
   proxy = @(X) proxZ(X,mu,opts);
   proxf = @(X) prox_tvz(X,proxz);
   proxg = @(X) prox_tvy(X,proxy);
   out=proxDykstra(Y,proxf,proxg,mit);
end

function X=proxZ(Z,lam,opts)
   out=runMethod(Z,lam,opts);
   X=out.X;
end
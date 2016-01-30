function out=proxDykstra(Y,proxf,proxg,mit)
% PROXDYKSTRA -- meta solvero
%
% min 0.5|X-Y|^2 + f(X) + g(X)
%
% In:
%   Y -- vector, matrix, tensor, whatever
%   proxf -- prox operator for f(X)
%   proxg -- proxg operator for g(X)
%
% Author: Suvrit Sra
%
% See also: DR

   X=Y;
   P=zeros(size(X));
   Q=zeros(size(Y));
   out.time0 = cputime;

   for i=1:mit
     Z=proxf(X+P);
     P=X+P-Z;
     Xold=X;
     X=proxg(Z+Q);
     Xp = permute(X, [2 3 1]);

     out.time(i) = cputime;
     out.delta(i) = norm(Xold(:)-X(:))/sqrt(prod(size(X)));
     out.snr(i) = snr(Xp);
     % TODO(mwytock): Remove hardcoded lambda
     out.objval(i) = 0.5*norm(X(:)-Y(:))^2 + 0.2*tv_image(Xp);
     Q=Z+Q-X;

     fprintf('%3d\t%4.2f\t%4.2f\n', i, out.snr(i), out.objval(i));
   end
   out.X = X;
end
function out=proxDR(Y,proxf,proxg,opt)
% PROXDR -- meta solver
%
% min 0.5|X-Y|^2 + f(X) + g(X)
%
% In: 
%   Y -- vector, matrix, tensor, whatever
%   proxf -- prox operator for f(X)
%   proxg -- proxg operator for g(X)
%
% 
% Author: Suvrit Sra
%  
% See also: proxDykstra
   
   X=Y;
   for i=1:opt.maxit
      Rg=2*proxg(X)-X;                   % Rg(X)
      X=proxf(Rg)-0.5*Rg + 0.5*X;        % 0.5*Rf(Rg(X))+0.5X
   end
   X=proxg(X);
end
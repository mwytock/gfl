function C = leftmultiplybyXt(Y,w)
% Fast computation for X'*Y for the group fused Lasso
%
% C = leftmultiplybyXt(Y,w)
%
% Compute X'*Y where X is the n*(n-1) design matrix for the weighted group
% fused Lasso, with weights defined by the vector w, and Y is any n*p
% matrix. The computation is done in O(np).
%
% Input
% Y :   a n*p matrix
% w :   (n-1)*1 vector of weights
%
% Output
% C :   the (n-1)*p matrix equal to X'*Y
%
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

% See paper for justification of the following formula
[n,p]=size(Y);
u=cumsum(Y);
C = ([1:n-1]'*u(end,:)/n - u(1:end-1,:)).*w(:,ones(1,p));

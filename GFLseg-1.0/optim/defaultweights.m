function w = defaultweights(n)
% Compute default weights for group fused Lasso
%
% w = defaultweights(n)
%
% The default weight is w(i)=sqrt(n/(i*(n-i)).
%
% Input
% n :   length of the signal
%
% Output
% w :   1*(n-1) vector of weights
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

a=[1:n-1]';
w = sqrt(n./(a.*(n-a)));

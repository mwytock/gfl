function g = XtX(n,A,B,w)
% Compute a submatrix of X'*X for the group fused Lasso
%
% g = XtX(n,A,B,w)
%
% Compute X(:,A)'*X(:,B), where X is the design matrix of the weigthed
% graph fused Lasso.
%
% X is a n*(n-1) matrix whose j-th column has the value:
% X(i,j) = w(j)*(j/n-1) for 1 <= j <= i
% X(i,j) = w(j)*j/n     for i+1 <= j <= n
%
% INPUT
% n :   size of the graph fused Lasso problem (X is a n*(n-1) matrix)
% A :   a column vector of indices in [1,n-1]
% B :   a column vector of indices in [1,n-1]
% w :   (n-1)*1 column vector of weights
%
% OUTPUT
% C = X(:,A)'*X(:,B), a length(A)*length(B) matrix
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


% Trivial case
if length(A)==0
    g=[];
    return
end

% Compute g(i,j) = min(A(i),B(j)) * (n - max(A(i),B(j))) * w(A(i)) *
% w(B(j)) / n 
u = A(:,ones(1,length(B)));
v = B(:,ones(1,length(A)))';
g = min(u,v).*(n-max(u,v)).*(w(A)*w(B)')/n;

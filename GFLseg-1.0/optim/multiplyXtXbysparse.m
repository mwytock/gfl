function r = multiplyXtXbysparse(n,ind,val,w)
% Fast computation of X'*X*b when b is sparse, for the group fused Lasso
%
% r = multiplyXtXbysparse(n,ind,val,w)
%
% Compute r = X'*X*b , where b is a row-sparse (n-1)*p matrix and X is the
% n*(n-1) design matrix for the weighted group fused lasso.
%
% INPUT
% n :   size of the problem
% ind : &*1 vector of indices of the non-zero rows of b (each in [1,n-1])
% val : a*p matrix whose rows are the non-zero rows of b (same order as
% ind)
% w :   (n-1)*1 vector of weights
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

[a p] = size(val);

if a==0
    % multiplication by zero: return 0
    r=zeros(n-1,p);
else
    
    % First multiply beta by the weights
    val = val.*w(ind,ones(1,p));
    
    % compute the matrix s of increments of r
    s=zeros(n-1,p);
    
    s(ind,:)=val;
    s=flipud(cumsum(flipud(s)));
    u=ind'*val;
    s=s-u(ones(n-1,1),:)/n;
    
    % then make the cumsum
    r=cumsum(s);

    % then multiply the rows by the weights
    r=r.*(w(:,ones(1,p)));
end

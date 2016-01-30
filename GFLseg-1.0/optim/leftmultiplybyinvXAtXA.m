function r = leftmultiplybyinvXAtXA(n,ind,val,w)
% Fast computation of inv(X'*X)*b for the group fused Lasso
%
% r = leftmultiplybyinvXAtXAweighted(ind,b,n,w)
%
% Compute r = inv(X(:,ind)'*X(:,ind))*b , where X is the n*(n-1) design
% matrix for the weighted group fused lasso. 
%
% INPUT
% ind : a*1 vector of indices between 1 and n-1, sorted in increasing order
% b :   a*p matrix
% n :   the size of X is n*(n-1)
% w :   (n-1)*1 vector of weights
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

[a p] = size(val);

if a==0
    % multiplication by zero: return 0
    r=zeros(a,p);
else
    % see paper for explanation of this formula
    u=diff([0;ind;n]);
    val=val./w(ind,ones(1,p));
    delta = diff([zeros(1,p);val;zeros(1,p)])./u(:,ones(1,p));
    r=-diff(delta);
    r=r./w(ind,ones(1,p));
end

function d=inftwonorm(X)
% ONETWONORM -- weighted 1,2 norm
   d=sum(X.^2,1);d=full(d);d=sqrt(d);
   d=max(d);
end
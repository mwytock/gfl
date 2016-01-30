function f = tv_image(X)
[m n k] = size(X);
assert(k==3);

f = 0;
for i=1:m
  f = f + sum(norms(X(i,2:end,:)-X(i,1:end-1,:),[],3));
end
for j=1:n
  f = f + sum(norms(X(2:end,j,:)-X(1:end-1,j,:),[],3));
end

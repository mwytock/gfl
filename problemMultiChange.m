
function Y = problemMultiChange(T,n,k)
% Generate problem instance with k segments
Y = zeros(n,T);
for i=1:k
  ii = ceil((i-1)/k*T)+1:ceil(i/k*T);
  Y(:,ii) = repmat(randn(n,1),1,length(ii))+0.1*randn(n,length(ii));
end

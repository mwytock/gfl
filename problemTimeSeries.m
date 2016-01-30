function prob = problemTimeSeries(T,p,a,e)
% Generate an AR(p) process with T time points and a switch in the middle

A1 = [randn(1,p)*a(1); eye(p-1) zeros(p-1,1)];
A2 = [randn(1,p)*a(2); eye(p-1) zeros(p-1,1)];

delta1 = max(real(eig(A1)));
delta2 = max(real(eig(A2)));
if delta1 > 1
  A1 = A1 - (delta1-1+1e-3)*eye(p);
end
if delta2 > 1
  A2 = A2 - (delta2-1+1e-3)*eye(p);
end

% General state space model
X = zeros(T+1,p);
X(1,:) = randn(p,1);
prob.A = sparse(T,T*p);
for t=2:T+1
  A = (t<T/2)*A1 + (t>=T/2)*A2;
  X(t,:) = A*X(t-1,:)' + e*randn(p,1);
  prob.A(t-1,(t-2)*p+1:(t-1)*p) = X(t-1,:);
end

prob.A1 = A1;
prob.A2 = A2;
prob.y = X(2:T+1,1);
prob.A = prob.A;
function [obj,gap] = primal_objval2(Y,U,lam)
V = projInfTwoBall(lam,U);
X = Y-rightMultByEtrans(V);
XE = X(:,1:end-1) - X(:,2:end);
gap = lam*onetwonorm(XE)-sum(sum((V.*XE)));
obj = 0.5*norm(X-Y,'fro')^2 + lam*onetwonorm(XE);
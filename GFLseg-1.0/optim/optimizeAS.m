function [beta AS out] = optimizeAS(beta0, AS0, lambda, XtY, Y, option)
% Solve the group fused Lasso optimization problem
%
% [beta AS] = optimizeAS(beta0, AS0, lambda, XtY, option)
%
% Optimize beta with an active set strategy. beta0 and AS0 are the initial
% guess for beta and the active set.
%
% INPUT
% AS0 :     a*1 vector with the indices of the blocks in the active set
% beta0 :   a*p matrix to initiliaze the optimization
% lambda :  the regularization parameter lambda
% XtY :     (n-1)*p matrix X'*Y
%
% OUTPUT
% beta :    b*p matrix of solutions
% AS :      indices of the active set, ie, of the non-zero rows of beta
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

% Initialization
AS=AS0;
beta=beta0;
[n p] = size(XtY);
n=n+1;
globalSol = 0;

t_start = cputime;
iter = 1;
meansignal = mean(Y);

% Main loop
while ~globalSol
    % Optimize over the current active set only
    beta = blockoptimizeAS( beta, AS, lambda, XtY(AS,:), n, option);

    % Remove from active set the zero coefficients
    nonzerocoef = sum(beta.^2,2)~=0;
    AS = AS(nonzerocoef);
    beta = beta(nonzerocoef,:);


    % Check optimality
    S =  multiplyXtXbysparse(n,AS,beta,option.weights)-XtY;
    normS = sum(S.^2,2);

    if isempty(AS)
        [maxnormS , imax] = max(normS);
        if maxnormS<lambda^2+option.tol
            % the optimal solution is the null matrix
            globalSol=1;
        else
            % this should only occur at the first iteration
            AS = imax;
            beta = zeros(1,p);
        end

    else

        % At optimality we must have normS(i)=lambda^2 for i in AS and
        % normS(i)<lambda^2 for i not in AS.

        lagr = max(normS(AS)); % This is equal to lambda^2 in theory. In practice, we use min(lambda,lagr).
        lagr = min(lagr,lambda^2);
        nonas = setdiff([1:n-1],AS);
        % Check optimality conditions for blocks not in the active set
        [y , i] = max(normS(nonas));
        if (isempty(nonas) || (y<lagr+option.tol))
            % Optimality conditions are fulfilled: we have found the global
            % solution
            globalSol = 1;
        else
          % Otherwise we add the block that violates most the optimality
          % condition
          AS = [AS ; nonas(i)];
          beta = [beta ; zeros(1,p)];
        end
    end


    out.iterTime(iter) = cputime;

    X = reconstructsignal(n,AS,beta,ones(n-1,1),meansignal);
    XE = X(2:end,:) - X(1:end-1,:);
    out.obj(iter) = 0.5*norm(X-Y,'fro')^2 + lambda*onetwonorm(XE');
    out.X = X;

    if mod(iter, 10)==0 && option.verbose >= 1
      fprintf('%d %f\n', iter, out.obj(iter));
    end
    out.statsTime(iter) = cputime;
    if cputime-t_start > option.max_time
      break
    end
    iter = iter + 1;
end
out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;

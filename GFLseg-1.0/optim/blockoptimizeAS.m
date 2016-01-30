function beta = blockoptimizeAS( beta0, AS, lambda, XtY, n, option)
% Block coordinate descente over the active set for group fused Lasso
%
% beta = blockoptimizeAS( beta0, AS, lambda, XtY, n, option)
%
% INPUT
% beta0 :  a*p matrix to initiliaze the optimization
% AS :     a*1 vector with the indices of the blocks in the active set
% lambda : the regularization parameter lambda
% XtY :    a*p matrix X'*Y(AS,:)
% n :      length of profiles
% option.tol : stop iterations when the gain on each block has been below
% this value [default=1e-8]
% option.maxit : maximum number of iterations [defautl=1e5]
% option.verbose : display information if 1
% option.weights : the weights of weighted group fused lasso
%
% OUTPUT
% beta :  the a*p matrix after optimization
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


% Print information
if option.verbose
    fprintf('\nEntering blockoptimize AS with %d active blocks\nmaxit=%d tol=%g verbose=%d\n',length(AS),option.maxit,option.tol,option.verbose);
end

% Prepare variables
beta = beta0;
[a,p] = size(beta);
normbeta = sqrt(sum(beta.^2,2));
tol = option.tol*p; % Note that we scale the tolerance tol by the dimension p
gain = 2*tol*ones(a,1);


% MAIN LOOP
if a>0
    % We optimize each block in turn
    itc=0;  % Iteration counter

    while(any(gain>tol) && (itc < option.maxit))

        % The loop continues as long as the gain in objective function
        % during the last optimization of at least one block has been above
        % the threshold tol, but we stop after a maximum of maxit
        % iterations.

        i = mod(itc,a)+1; % index of block to update in AS
        ASi = AS(i);      % block to update

        % Compute the vector S
        XitX = XtX(n,ASi,AS,option.weights);  % compute XitX = X(:,ASi)'*X(:,AS)
        gammai = XitX(i);      % compute gammai = X(:,i)'*X(:,i) (= ASi*(n-ASi)/n)
        indwithouti = [1:i-1,i+1:a]; % the indices of the active set excluding the one being optimized
        S = XtY(i,:) - XitX(indwithouti)*beta(indwithouti,:); % compute S=X(:,ASi)'*Y - X(:,ASi)*X(:,AS(indwithouti))*beta(AS(indwithouti),:)

        % Check the norm of S to decide where the optimum is
        nS = norm(S);
        if (nS < lambda)
            newbeta = zeros(1,p); % If ||S||<lambda then the optimal beta(ASi,:) is zero
        else
            newbeta = S*(1-lambda/nS)/gammai; % Otherwise the optimal beta(ASi,:) is obtained by shrinking S
        end

        % Compute the gain in the objective function at this iteration
        newnormbeta = norm(newbeta);
        gain(i) = (gammai*(normbeta(i) + newnormbeta)/2 + lambda)*(normbeta(i)-newnormbeta) + S*(newbeta' - beta(i,:)');
        if option.verbose >= 2
            fprintf('[blockoptimizeAS] Iteration %d, update block %d, gain=%g\n',itc,ASi,gain(i));
        end

        % Update beta
        beta(i,:) = newbeta;
        normbeta(i) = newnormbeta;

        % increase loop iteration counter
        itc=itc+1;
    end
end

function res = gflars(Y,k,option)
% Segmentation of a multi-dimensional signal with the group fused LARS. 
%
% res = gflars(Y,k)
% res = gflars(Y,k,option)
%
% INPUT
% Y :       a n*p signal to be segmented
% k :       the number of change points to find
% option :  an optional list of parameters:
%   - option.epsilon : values smaller than epsilon are considered null
%   [default=*1e-9] 
%   - option.verbose : display information if 1
%   - option.weights : a (n-1)*1 vector of weights for the weigthed graph
%   fused Lasso penalty. If absent, the default weights sqrt(n/(i*(n-i)))
%   are taken.
%
% OUTPUT
% res.lambda : the estimated lambda values for each change-point
% res.jump :   the successive change-point positions (1*k)
% res.value{i} :  a i*p matrix of change-point values for the first i
% change-points
% res.meansignal: the mean signal per column (1*p vector)
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


if nargin<2
    fprintf('Error: too few arguments\nUsage:\nres = gflars( Y , k [,option] )\n');
    return;
end


% Size of the signal
% n is the length of the signal
% p is the dimension of the signal (or number of profiles)
[n p]=size(Y);

% Read options

if nargin==2
    option=[];
end

if ~isfield(option,'epsilon')
    option.epsilon = 1e-9;
end

if ~isfield(option,'verbose')
    option.verbose = 0;
end

if isfield(option,'weights')
    weights = option.weights;
else
    weights = defaultweights(n);
end

if option.verbose
    fprintf('\nStarting gflars with epsilon=%g verbose=%d\n',option.epsilon,option.verbose);
end







% Store the mean signal per column (useful for reconstruction)
res.meansignal = mean(Y);

res.lambda=zeros(k,1);
res.jump=zeros(k,1);
EPSILON = option.epsilon;

% Initialize cHat = X'*Y
cHat = leftmultiplybyXt(Y,weights);

% Main loop to find the successive jumps
for iter=1:k
        
    % Compute the row norms of cHat
    cHatSquareNorm = sum(cHat.^2,2);
    [bigcHat,besti]=max(cHatSquareNorm);

    % In the first iteration, we add the most correlated feature to the
    % active set. For the other iterations, this is already done at the end
    % of the previous iteration
    if iter==1
        res.jump(1)=besti;
    end
    
    % Compute the descent direction w = inv(X(:,A)'*X(:,A))*cHat(A,:)
        
    [A,I]=sort(res.jump(1:iter)); % TODO: could be done faster. Only required to call the function leftmultiplybyinvXAtXA
    w = leftmultiplybyinvXAtXA(n,A,cHat(A,:),weights);
    
    % Comute a = X'*X*w
    a = multiplyXtXbysparse(n,A,w,weights);
    
    % Compute the descent step
       % For each i we find the largest possible step alpha by solving:
    % norm(cHat(i,:)-alpha*a(i,:)) = norm(cHat(j,:)-alpha*a(j,:))
    % where j is in the active set.
    % We write it as a second order polynomial 
    % a1(i)*alpha^2 - 2* a2(i)*alpha + a3(i)
    
    a1 = bigcHat - sum(a.^2,2);
    a2 = bigcHat - sum(a.*cHat,2);
    a3 = bigcHat - cHatSquareNorm;

    % We solve it
    gammaTemp = zeros(2*(n-1),1);
    
    % First those where we really have a second-order polynomial
    subset = find(a1 > EPSILON);
    gammaTemp(subset) = (a2(subset) + sqrt(a2(subset).^2 - a1(subset).*a3(subset)))./a1(subset);
    gammaTemp(subset+n-1) = (a2(subset) - sqrt(a2(subset).^2 - a1(subset).*a3(subset)))./a1(subset);
    
    % then those where the quadratic term vanishes and we have a
    % first-order polynomial
    subset = find((a1 <= EPSILON) & (a2 > EPSILON));
    gammaTemp(subset) = a3(subset) ./ (2*a2(subset));
    gammaTemp(subset+n-1) = a3(subset) ./ (2*a2(subset));

    % Finally the active set should not be taken into account, as well as
    % those for which the computation gives dummy solutions
    maxg=max(gammaTemp)+1;
    subset = find((a1 <= EPSILON) & (a2 <= EPSILON));
    gammaTemp(subset) = maxg;
    gammaTemp(n+subset) = maxg;
    gammaTemp(A) = maxg;
    gammaTemp(n+A-1) = maxg;    
    gammaTemp(gammaTemp<=0)=maxg;
    gammaTemp(imag(gammaTemp)~=0) = maxg;

    % Now we can take the minimum
    [gamma,nexttoadd]=min(gammaTemp);
    
    % Update
    res.value{iter} = zeros(iter,p);
    res.value{iter}(I,:) = gamma*w;
    if iter>1
        res.value{iter}(1:(iter-1),:) = res.value{iter}(1:(iter-1),:) + res.value{iter-1};
    end
    res.lambda(iter)=sqrt(bigcHat);
    
    if iter<k
        res.jump(iter+1) = 1+mod(nexttoadd-1,n-1);
        cHat = cHat-gamma*a;
    end
end


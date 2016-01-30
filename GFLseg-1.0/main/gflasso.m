function res = gflasso(Y,lambda,option)
% Segmentation of a multi-dimensional signal with the group fused Lasso.
%
% res = gflasso(Y,lambda)
% res = gflasso(Y,lambda,option)
%
% INPUT
% Y :       a n*p signal to be segmented
% lambda :  a vector of lambda values for the regularization parameter of
% the group Lasso.
% option :  an optional list of parameters:
%   - option.tol : stopping criterion on the objective function decrease
%   in the block coordinate descent [default=p*1e-8]
%   - option.maxit : maximum number of iterations in the block coordinate
%   descent [defautl=1e5]
%   - option.verbose : display information if 1
%   - option.weights : a (n-1)*1 vector of weights for the weigthed graph
%   fused Lasso penalty. If absent, the default weights sqrt(n/(i*(n-i)))
%   are taken.
%
% OUTPUT
% res.lambda :    the lambda values in decreasing order
% res.jump{i} :   a j*1 vector of change-point positions for the i-th
% lambda value (j depends on lambda)
% res.value{i} :  a j*p matrix of change-point values for the i-th lambda
% value
% res.meansignal: the mean signal per column (1*p vector)
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


if nargin<2
    fprintf('Error: too few arguments\nUsage:\nres = gflasso( Y , lambda [,option] )\n');
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

if ~isfield(option,'maxit')
    option.maxit = 1e5;
end

if ~isfield(option,'tol')
    option.tol = 1e-8;
end

if ~isfield(option,'verbose')
    option.verbose = 0;
end

if ~isfield(option,'weights')
    option.weights = defaultweights(n);
end


if option.verbose
    fprintf('\nStarting gflasso with maxit=%d tol=%g verbose=%d\n',option.maxit,option.tol,option.verbose);
end


% Store the mean signal per column (useful for reconstruction)
res.meansignal = mean(Y);

% Sort the lambda by decreasing value, in order to use the group Lasso path
% algorithm with active set.
lambda = sort(lambda,'descend');
nlambda = length(lambda);
res.lambda = lambda;


% Compute C = X'*Y
C = leftmultiplybyXt(Y,option.weights);

% % Prepare the result variable
% res.loss = zeros(nlambda,1);
% res.penalty = zeros(nlambda,1);
% normY2 = norm(Y)^2;

% Initialize the active set and the solution
as=[];
beta=zeros(0,p);


% Iterate over lambda
for i=1:nlambda

    % Optimize for this lambda with a warm restart from the previous lambda
    [beta as out] = optimizeAS(beta, as, lambda(i), C, Y, option);
    res.jump{i} = as;
    res.value{i} = beta;
    res.time = out.time;
    res.obj = out.obj;
    res.X = out.X';

%     % Compute losses etc
%     res.penalty(i) = sum(sqrt(sum(beta.^2,2))); % L1,2 norm of beta
%     res.loss(i) = normY2/2;
%     nas=length(as);
%     if nas>0
%         tmp = as*ones(1,nas);
%         xtxas = min(tmp,tmp')*(n-max(tmp,tmp'))/n; % X(:,AS)'*X(:,AS)
%         res.loss(i) = res.loss(i) + sum(sum(beta .* (xtxas*beta)))/2 - sum(sum(C(as,:).*beta));
%     end

end

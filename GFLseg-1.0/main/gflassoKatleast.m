function res = gflassoKatleast(Y,k,option)
% Group fused Lasso segmentation with at least k change-points
%
% res = gflassoKatleast(Y,k)
% res = gflassoKatleast(Y,k,option)
%
% Segmentation of a multi-dimensional signal with the group fused Lasso.
% Runs along the path, with decreasing lambda, until at least k
% change-points are found. lambda starts with the maximum value which leads
% to at least one change-points, and is then iteratively divided by a decay
% factor.

% INPUT
% Y :   n*p signal to be segmented
% k :   number of jumps  
% option :  an optional list of parameters:
%   - option.decay : the rate of geometric decay of lambda on the path
%   (default:0.5)
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
% res.lambda : regularization paramter for k jumps
% res.jump : a k*1 vector of change-point positions
% res.value : a k*p matrix of change-point values 
% res.meansignal: the mean signal per column (1*p vector)
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert



if nargin<2
    fprintf('Error: too few arguments\nUsage:\nres = gflassoK( Y , k [,option] )\n');
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

if ~isfield(option,'decay')
    option.decay = 0.5;
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
    fprintf('\nStarting gflassoK with decay=%g maxit=%d tol=%g verbose=%d\n',option.decay, option.maxit,option.tol,option.verbose);
end



% Store the mean signal per column (useful for reconstruction)
res.meansignal = mean(Y);

% Compute C = X'*Y
C = leftmultiplybyXt(Y,option.weights);

% The maximum lambda value (that leads to a jump)
lambdamax = max(sqrt(sum(C.^2,2)))*(1-option.tol);

% Initialize the active set and the solution
as=[];
beta=zeros(0,p);

% We can not find more than n-1 change-points
k=min(k,n-1);

lambda = lambdamax/option.decay; % just a trick to start at lambdamax later

while length(as)<k

    % Decrease lambda
    lambda = lambda*option.decay;
    
    % Optimize for the current lambda with a warm restart from the previous
    % solution
    
    if option.verbose
        fprintf('Now trying lambda=%g ...',lambda);
    end
    
    [beta as] = optimizeAS(beta, as, lambda, C, option);
    
    if option.verbose
        fprintf('Found %d jumps\n',length(as));
    end
    
end

res.lambda = lambda;
res.jump = as;
res.value = beta;

function res = gflassoK(Y,k,option)
% Group fused Lasso segmentation of a multi-dimensional signal with k change-points 
%
% res = gflassoK(Y,k)
% res = gflassoK(Y,k,option)
%
% This script automatically finds a regularization parameter that leads to
% a given number of jumps, by dichotomic search.
%
% INPUT
% Y :   n*p signal to be segmented
% k :   number of jumps  
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
% OUTPIT
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
    fprintf('\nStarting gflassoK with maxit=%d tol=%g verbose=%d\n',option.maxit,option.tol,option.verbose);
end



% Store the mean signal per column (useful for reconstruction)
res.meansignal = mean(Y);

% Compute C = X'*Y
C = leftmultiplybyXt(Y,option.weights);

% The maximum lambda value (that leads to 0 jump)
lambdamax = max(sqrt(sum(C.^2,2)))*(1+option.tol);
lambdamin = 0;
MINLAMBDADIFF = lambdamax/1e10;

% Initialize the active set and the solution
as=[];
beta=zeros(0,p);

% Dichotomic search
while (length(as)~=k) && (lambdamax-lambdamin>MINLAMBDADIFF)
    
    lambda = (lambdamin+lambdamax)/2;
    
    % Optimize for this lambda with a warm restart from the solution at
    % lambdamax
    
    if option.verbose
        fprintf('Now trying lambda=%g lambdamax=%g lambdamin=%g\n',lambda,lambdamax,lambdamin);
    end
    
    [betanew asnew] = optimizeAS(beta, as, lambda, C, option);
    
    if option.verbose
        fprintf('Found jumps=%d\n',length(asnew));
    end
    
    if length(asnew) <= k
        beta = betanew;
        as = asnew;
        lambdamax = lambda;
    else
        lambdamin = lambda;
    end
    
end

res.lambda = lambda;
res.jump = as;
res.value = beta;

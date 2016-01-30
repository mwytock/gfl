function res = dpseg(Y,option)
% Segmentation of a multi-dimensional signal with dynamic programming. 
% 
% res = dpseg(Y)
% res = dpseg(Y,option)
%
% INPUT
% Y :       a n*p signal to be segmented
% option :  a list of parameters:
%   - option.candidatechangepoints : a vector of candidate positions for
%   change-points (default=[1:n-1])
%   - option.kmax : maximum number of change-points to test (default :
%   length(option.candidatechangepoints)
%   - option.threshold : stopping criteria. Typically chosen to be in the interval 
% (0 0.5]. The smaller the threshold, the higher the tendency to keep more
% breakpoints. The criteria is based on the method found in 'Picard et al
% (2005)', "A statistical approach for array CGH data analysis" (BMC
% Bioinformatics). Default=0.5
%
% OUTPIT
% res.jump{i} :   a j*1 vector of change-point positions for the i-th
% lambda value (j depends on lambda). i varies between 1 and kmax
% res.rse : a (kmax+1)-dimensional vector of residual squared error
% res.kbest: the number of selected change-points
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert




% Size of the signal
% n is the length of the signal
% p is the dimension of the signal (or number of profiles)
[n p]=size(Y);

% Read options
if nargin==1
    option=[];
end

if ~isfield(option,'candidatechangepoints')
    option.candidatechangepoints = [1:n-1];
end

kmaxmin = floor(n/4); % maximum number of change-points we can detect

if ~isfield(option,'kmax')
    option.kmax = min(length(option.candidatechangepoints),kmaxmin);
end

% Do not try to find more than kmaxmine change-points (otherwise the likelihood
% may diverge etc...)
if option.kmax > kmaxmin
    fprintf('Warning: not enough points to optimize the number of change-points up to %d\n',option.kmax)
    option.kmax = kmaxmin;
    fprintf('Set the maximum number of change-points to %d\n',option.kmax)
end

if ~isfield(option,'threshold')
    option.threshold = 0.5;
end

% Compute boundaries of the smallest intervals considered. 
b = sort(union([0],union([n],option.candidatechangepoints)));
k = length(b)-1; % k is the number of such intervals

% Compute the k*k matrix J such that J(i,j) for i<=j is the RSE when
% intervales i to j are merged
J = zeros(k,k);
s=[zeros(1,size(Y,2));cumsum(Y)]; % cumsum of the rows
v=[0;cumsum(sum(Y.*Y,2))]; % cumsum of squared norm of the rows

for i=1:k
    for j=i:k
        Istart=b(i)+1; 
        Iend=b(j+1); 
        J(i,j) = v(Iend+1)-v(Istart) - sum((s(Iend+1,:)-s(Istart,:)).^2)/(Iend-Istart+1);
    end
end

% Dynamic programming recursion
V=zeros(option.kmax+1,k); % V(i,j) is the best RSE for segmenting intervals 1 to j with at most i-1 change-points
jump=zeros(option.kmax,k);
% With no change-points, V(1,j) is just the precomputed RSE for intervals 1
% to j
V(1,:)=J(1,:);

% Then we apply the recursive formula
for ki=1:(option.kmax)
    for j=(ki+1):k
        [val,ind] = min(V(ki,ki:j-1) + J(ki+1:j,j)');
        V(ki+1,j) = val;
        jump(ki,j) = ind+ki-1;
    end
end

% Optimal segmentations
for ki=1:option.kmax
    res.jump{ki} = zeros(1,ki);
    res.jump{ki}(ki) = jump(ki,k);
    for i=ki-1:-1:1
        res.jump{ki}(i) = jump(i,res.jump{ki}(i+1));
    end
end

% Convert back the index of the interval to the last position before the
% jump
rightlimit = b(2:end);
for ki=1:option.kmax
    res.jump{ki} = rightlimit(res.jump{ki});
end

% RSE as a function of number of change-points
res.rse=V(:,k);

% Optimal number of changepoints
J = log(res.rse); % log-likelihood is -n/2*(J-log(n)+1+log(2*pi));
Km=length(J);
Jtild=(J(Km)-J)/(J(Km)-J(1))*(Km-1)+1; %normalize
res.kbest=max(find(diff(diff(Jtild))>option.threshold))+1; % find inflexion point
if isempty(res.kbest)==1 res.kbest=1; end;

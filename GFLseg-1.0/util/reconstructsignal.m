function signal = reconstructsignal(n,ind,val,w,meansignal)
% Reconstruct a piecewise-constant profile from its increments.
%
% p = reconstructsignal(n,ind,val,w)
% p = reconstructsignal(n,ind,val,w,meansignal)
%
%
% INPUT
% ind :     a*1 vector of change-point positions. A position i means that the
% change-point occurs between i and i+1. Positions must therefore be in the
% range [1,n-1], where n is the length of the signal
% val :     a*p matrix of change-point values. The i-th row is the
% increment of the i-th change-point at position jump(i) (which .
% n :       length of the signal
% w :       (n-1)*1 vector of weights (in the weighted group fused Lasso)
% meansignal (optional): 1*p vector of mean signal across columns (default
% zeros)
%
% OUTPUT
% signal : the n*p signal with corresponding change-points, and such that
% mean(signal)=meansignal.
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


[a,p]=size(val);
% a : number of change-points
% p : dimension of the signal

% If meansignal is not provided as argument, set it to zero
if nargin<5
    meansignal=zeros(1,p);
end

% Reconstruct the signal
signal = zeros(n,p);
signal(ind+1,:) = val.*w(ind,ones(1,p));
signal = cumsum(signal);

% Set the mean signal to meansignal
m = meansignal - mean(signal);
signal = signal + m(ones(n,1),:);

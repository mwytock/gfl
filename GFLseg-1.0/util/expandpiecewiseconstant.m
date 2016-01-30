function Y = expandpiecewiseconstant(jumps,val)
% Expand a piecewise constant signal from a compact to an explicit form
%
% Y = expandpiecewiseconstant(jumps,val)
%
% We usually store a p-dimensional piecewise constant signal of length n
% with k constant intervals in a compact form, by only storing the last
% position of each interval and the p-dimensional value taken on each
% interval. This function expands this compact representation in an
% explicit n*p vector.
%
% INPUT
% jumps :  the k positions of ends of the intervals. By convention, the
% last position n=jumps(end) is the end of the signal
% val :    a k*p matrix are the values on each interval
%
% OUTPUT
% Y :      a n*p matrix of piecewise constant signals
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

% Number of intervals and dimension
[k,p] = size(val);

% Boundaries of the intervals ( between b(i)+1 and b(i+1) )
b = [0;jumps];

Y = zeros(b(end),p);

% Expand Y on each interval
for i=1:k
    Istart = b(i)+1;
    Iend = b(i+1);
    Y(Istart:Iend,:) = val(i*ones(Iend-Istart+1,1),:);
end
    
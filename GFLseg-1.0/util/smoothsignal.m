function s = smoothsignal(Y,jumps)
% Create a piecewise constant profile by averaging Y between jumps.
%
% s = smoothsignal(Y,jumps)
%
% INPUT
% Y : a n*p signal (length n, dimension p)
% jumps : a vector of jump positions. j(i) is the last position of the i-th
% interval.
%
% OUTPUT
% s : a structure with s.jumps the list of k jumps (including n) and s.val
% the k*p matrix of values over the k intervals.
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

% Size of the signal
[n,p]=size(Y);

% Boundaries of the intervals (i-th interval from b(i)+1 to b(i+1) )
b=sort(union([0],[jumps;n]));

% Number of intervals
k = length(b)-1;

% Jumps (including n) defined as the last position of each interval
s.jumps = b(2:end);

% Values on each interval
s.val = zeros(k,p);
for i=1:k
    Istart = b(i)+1;
    Iend = b(i+1);
    s.val(i,:) = mean(Y(Istart:Iend,:));
end
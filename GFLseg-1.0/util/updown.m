function ud = updown(Y)
% Detect frequent gains an losses
%
% ud = updown(Y)
%
% Compute statistics to assess significant positive and negative values on
% each row of Y. Currently, compute the mean of positive values for the
% "up" statistics, the mean of the negative values for the "down"
% statistics.
%
% Input:
% Y : a n*p signal (length n, dimension p)
%
% Output:
% ud : a n*2 matrix, the first column are the "up" statistics, the second
% are the "down" statistics
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

[n,p] = size(Y);

ud = zeros(n,2);
for i=1:n
    vali = Y(i,:);
    ud(i,1) = mean(vali(vali>0));
    ud(i,2) = mean(vali(vali<0));
end
ud(isnan(ud))=0;

    

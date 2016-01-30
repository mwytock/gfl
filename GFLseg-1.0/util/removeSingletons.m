function r = removeSingletons(res)
% Remove intervals of length 1 from a piecewise constant signal.
% 
% INPUT
% res : the output of simpleGFL, which codes for a piecewise constant
% function with the fields res.jumps (last position of each interval),
% res.val (values on each interval), res.updown (up/down score on each
% interval)
%
% OUTPUT
% r : same as the input, with intervals of lenght 1 deleted
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

% Length of each interval
d = [res.jumps(1) ; diff(res.jumps)];

% Keep only the ones longer than 1
keepme = d>1;

r.jumps = res.jumps(keepme);
r.val = res.val(keepme,:);
r.updown = res.updown(keepme,:);

% In case the last intervals are removed, we extend the last interval kept
% to the total length
r.jumps(end) = res.jumps(end);
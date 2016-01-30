function y=randomprofile(LENGTH,NJUMPS,NOISELEVEL,DIM)
% Generate a random DIM-dimensional profile with jumps and noise
%
% y=randomprofile(LENGTH,NJUMPS,NOISELEVEL,DIM)
%
% Generate a random profile (vector) of length length, with NJUMPS jumps
% randomly chosen. Between two jumps, the profile is constant, uniformly
% chosen between 0 and 1, and a Gaussian noice of variance NOISELEVEL is
% added.
%
% y is a structure that contains two fields:
%   y.profile is the profile (a LENGTH by DIM vector)
%   y.jumps is the list of jump positions (the last position on the left of
%   a jump)
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

% First make the noise
y.profile=randn(LENGTH,DIM)*NOISELEVEL;

% Choose the jumps
u=randperm(LENGTH-1);
y.jumps = sort(u(1:NJUMPS));
jumps=[0 y.jumps LENGTH];

% Add the random piecewise linear profile
for i=1:(NJUMPS+1)
    y.profile(jumps(i)+1:jumps(i+1),1:DIM) = y.profile(jumps(i)+1:jumps(i+1),1:DIM) + ones(jumps(i+1)-jumps(i),1)*rand(1,DIM);
end

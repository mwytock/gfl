function res = simpleGFL(Y,option)
% Automatic multiple change-point detection
%
% res = simpleGFL(Y,option)
%
% This function is a wrapper for signal segmentation into at least k
% change-points by the group fused Lasso or LARS, followed by change-point
% selection with a dynamic programming optimization.
%
% INPUT
% Y :  a matrix to segment (each column is a signal)
% option.GFLalgo: one of 'GFlars' or 'GFlasso', the algorithm used for the
% first segmentation (default: 'GFlars')
% option.k :  the number of change-points to find in the first segmentation
% with group fused Lars or Lasso
% option.threshold : thresholf for dynamic programming (see DPSEG)
%
% OUTPUT
% res.jumps : the positiongs of the change-points 
% res.val : the value on each interval
% res.updown : the up and down statistics on each interval
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

if nargin == 1
    option=[];
end

if ~isfield(option,'GFLalgo')
    option.GFLalgo = 'GFlars';
end

if ~isfield(option,'k')
    option.k = min(100,size(Y,1)-1);
end



% Run GFLasso or GFLars

if strcmp(option.GFLalgo,'GFlars') == 1
    fprintf('Running gflars...');
    tic
    res1 = gflars(Y,option.k,option);
    toc
    
elseif strcmp(option.GFLalgo,'GFlasso') == 1
    fprintf('Running gflasso\n');
    res1 = gflassoKatleast(Y,option.k,option);
    
end



% Post-processing by dynamic programming
option.candidatechangepoints = res1.jump;
res2 = dpseg(Y,option);

% Extract piecewise constant approximation
res = smoothsignal(Y,res2.jump{res2.kbest});

% Compute up/down statistics
res.updown = updown(res.val);

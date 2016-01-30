function r = chrInfo(chrIndex)
% Compute several informations about chromosomes
%
% r = chrInfo(chrIndex)
%
% INPUT
% chrIndex : a vector of chromosome numbers (associated to each probe)
%
% OUTPUT
% r.chrs :         sorted list of chromosomes
% r.nchrs :        number of chromosomes
% r.chr_nums :     indices of the end of each chromosome
% r.chr_data_len : number of probes in each chromosome
% r.probes :       probes of each chromosome
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

r.chrs = unique(chrIndex); % sorted list of chromosomes
r.nchrs = length(r.chrs); % number of chromosome
r.chr_nums = zeros(1, r.nchrs); % indices of the end of each chromosomes
r.chr_data_len = zeros(1, r.nchrs); % length of each chromosomes
r.probes = cell(r.nchrs,1); % the probes of each chromosome

for c=1:r.nchrs
    % Identify the probes in this chromosome
    tmp = find(chrIndex==r.chrs(c));
    % Find the index of the last probe
    r.chr_nums(c) = tmp(end);
    % Number of probes
    r.chr_data_len(c) = length(tmp);
    % Store the probes
    r.probes{c} = tmp;
end
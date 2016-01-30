function [X,Y,chrIndex,clones,starts,stops] = importDataCGH(clonefile,datafile,I)
% Read multiple samples or array CGH data
%
% [X,Y,chrIndex,clones,starts,stops] = importDataCGH(clonefile,datafile)
% [X,Y,chrIndex,clones,starts,stops] = importDataCGH(clonefile,datafile,I)
%
% INPUT
% clonefile : a single column text file containing the clone names that
% correspond to the rows in the datafile
% datafile : a tab delimited text file containing informations about
% multiple CGH arrays. It should contain the following columns
%   1. Chromosome number
%   2. Clone start position (bp)
%   3. Clone end position (bp)
%   4. Log ratio sample 1
%   5. Log ratio sample 2
%   ...
% I: the indices of the samples we want to extract (default: all)
%
% OUTPUT
% X : position of each probe (sorted by increasing chromosome, and
% increasing position within each chromosome)
% Y : log-ratio of each probe on each sample
% chrIndex : the chromosome of each probe
% starts : start position of each probe
% stops :  stop position of each probe
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


% Read files
data = dlmread(datafile);
clones = textread(clonefile,'%s');

% By default extract all samples
if (nargin < 3)
    I = 1:(size(data,2)-3);
end

% We reorder the probes by increasing chromosome, and by increasing
% position within a chromosome
chr = data(:,1); % chromosome of each probe
chrs = unique(chr); % sorted list of chromosomes
nchrs = length(chrs); % number of chromosomes
pos = (data(:,2) + data(:,3) + 1)/2; % position of each probe

[n,p] = size(data);
sortind = zeros(n,1);
i=1;
% Loop over chromosomes (in increasing order)
for c = 1:nchrs
    ind = find(chr==chrs(c)); % probes on this chromosome
    nind = length(ind); % number of probes on this chromosome
    [u J] = sort(pos(ind));
    sortind(i:i+nind-1) = ind(J); % indices of the probes sorted by position
    i = i+nind;
end

% then we extract the informations we want
data = data(sortind,:); % sort probes by increasing chromosomes, and increasing position within each chromosome
clones = clones(sortind); % probe names
chrIndex = data(:,1); % chromosome of each probe
starts = data(:,2); % start position of the probes
stops  = data(:,3); % stop position of the probes
X = pos(sortind); % average position of the probes
Y = data(:,I+3); % matrix of log-ratio for samples in I

% Remove NaN
Y(isnan(Y)) = 0;

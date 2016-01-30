function gainlossCGH(clonefile,datafile,outfile,outpic,I)
% Detect gain and losses in CGH profiles
%
% gainlossCGH(clonefile,datafile,outfile)
% gainlossCGH(clonefile,datafile,outfile,I)
%
% INPUT
% clonefile : file that contains the names of the probes
% datafile : a tab delimited text file with the CGH profiles. See
% importData for the format of datafile
% outfile : a file name to print the result 
% I : a list of samples to consider in the data file (default: 1:n)
%
% OUTPUT
% This function produces two outputs:
%  1. a file 'outfile" in seeGH format that contains the 'up' and 'down'
%  statistics along the genome 
%  2. a figure showing the 'up' and 'down' statistics along the genome
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


% Read data
fprintf('## Reading the data\n')
if nargin<4
    [X,Y,chrIndex,clones,starts,stops] = importData(clonefile,datafile);
else
    [X,Y,chrIndex,clones,starts,stops] = importData(clonefile,datafile,I);
end

[n,p] = size(Y);
fprintf('Found %d profiles of length %d\n\n',p,n)

% Process chromosome information
chr = chrInfo(chrIndex);

% Run segmentation on each chromosome
fprintf('## Performing joint segmentation\n')
updown = zeros(n,2); % to store the up and down statistics
for c=1:chr.nchrs
    fprintf('Now segmenting chromosome %d..',c)
    res = simpleGFL(Y(chr.probes{c},:));
    % remove singletons (intervals of length 1)
    res = removeSingletons(res);
    updown(chr.probes{c},:) = expandpiecewiseconstant(res.jumps,res.updown);
end

% Print file with 'up' and 'down' statistics per probe
fprintf('\n## Print results in %s\n\n',outfile)
fid=fopen(outfile,'wt');
if fid == -1
    disp(['ERROR: Could not open ', outfile, ' for writing']);
else
    CLONES = strvcat(clones);
    for i=1:n
        fprintf(fid,'%s\t%d\t%d\t%d\t%1.4f\t%1.4f\n', ...
            CLONES(i,:),chrIndex(i), starts(i), stops(i), ...
            updown(i,1),updown(i,2));
    end
    fclose(fid);
end

% Pretty-plot the 'up-down' statistics
fprintf('## Plot a figure in %s\n',outpic)
figure;
plotUpDown(chr,updown,{'-g','-r'},[-1.5;1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script runs GFLseg on the bladder cancer dataset (see data directory for more
% information about the dataset).
% We perform joint automated segmentation of each chromosome, and detect
% regions of gain and loss from the smoothed profiles.
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Now running demoGFLlung.m')
disp('Automatic analysis of a set of bladder cancer CGH profiles')


%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%

% Set the directory containing the files
% (here we set '.' because the data files are already on the MALAB path)
% datadir = fullfile('..','data'); % directory that contain the files
datadir = '';

% Set the names of the two input files (clone names and log-ratio file)
% To see the format of these files, type help importDataCGH
clonefile = fullfile(datadir,'bladderClones.txt'); % names of each probe
datafile = fullfile(datadir,'bladderData.txt'); % data

% Load data
[X,Y,chrIndex,clones,starts,stops] = importDataCGH(clonefile,datafile);
[n,p] = size(Y);
fprintf('Found %d profiles of length %d\n\n',p,n)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform joint segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process chromosome information
chr = chrInfo(chrIndex);

% Run segmentation on each chromosome
updown = zeros(n,2); % to store the up and down statistics
smooth = zeros(n,p); % to store the smoothed profiles
njumps = 0; % total number of change-points

for c=1:chr.nchrs
    fprintf('Now segmenting chromosome %d..',c)
    % run segmentation
    res = simpleGFL(Y(chr.probes{c},:));
    njumps = njumps + length(res.jumps);
    % remove singletons (intervals of length 1)
%    res = removeSingletons(res);
    % reconstruct up/down statistics for this chromosome
    updown(chr.probes{c},:) = expandpiecewiseconstant(res.jumps,res.updown);
    % reconstruct smooth profiles for this chromosome
    s = smoothsignal(Y(chr.probes{c},:),res.jumps);
    smooth(chr.probes{c},:) = expandpiecewiseconstant(s.jumps,s.val);
end

fprintf('Total: %d change-points on %d chromosomes\n',njumps,chr.nchrs)

%%%%%%%%%%%%%
% Plot figure
%%%%%%%%%%%%%

figure;
% Plot a raw and smoothed profile
subplot(3,1,1)
plotChromosome(chr,Y(:,end),{'.b'})
hold on
plot(smooth(:,end),'k-','LineWidth',2)
hold off
% Plot all smoothed profiles
subplot(3,1,2)
plotChromosome(chr,smooth)
% Plot up/down statistics
subplot(3,1,3)
plotChromosome(chr,updown,{'-g','-r'},[-0.5;0.5])

disp('Look at the plot!')

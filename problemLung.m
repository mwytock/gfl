function Y = problemLung()
datadir='';

% Set the names of the two input files (clone names and log-ratio file)
% To see the format of these files, type help importDataCGH
clonefile = fullfile(datadir,'lungClones.txt'); % names of each probe
datafile = fullfile(datadir,'lungData.txt'); % data
I = 1:18; % Samples we want to extract

% Load data
[X,Y,chrIndex,clones,starts,stops] = importDataCGH(clonefile,datafile,I);
Y = Y';
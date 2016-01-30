function Y = problemBladder()

datadir = '';
clonefile = fullfile(datadir,'bladderClones.txt'); % names of each probe
datafile = fullfile(datadir,'bladderData.txt'); % data
% Load data
[X,Y,chrIndex,clones,starts,stops] = importDataCGH(clonefile,datafile);
Y = Y';

% This script compares the group fused Lasso and the fused Lasso
% on simulated data.
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

disp('Now running demoGFLtoy.m')
disp('A small change point detection problem with the fused Lasso and the group fused Lasso.')

% Initialize the random number generator (to be able to reproduce the same
% experiment)
s = RandStream.create('mt19937ar','seed',5490);
RandStream.setDefaultStream(s);

% Parameters of the experiments
n=500;          % length of profiles
p=3;            % number of profiles
k=5;            % number of change points
noiselevel=0.2;% noise level
nplot = min(p,3); % number of profiles to plot simultaneously

% Generate a random multi-dimensional profile with change points
y=randomprofile(n,k,noiselevel,p);

% Set weights for the (group) fused lasso
o.weights = defaultweights(n);


% Find change-points with group fused Lasso
disp('Now running gflasso')
tic;
reslasso = gflassoK(y.profile , k , o);
toc
% Plot results
figure(1)
clf
for ipro=1:nplot
    subplot(nplot,1,ipro);
    plot(y.profile(:,ipro),'o')
    hold on
    ylim=get(gca,'ylim');
    jumps = reslasso.jump';
    line([jumps;jumps],ylim,'linewidth',2,'color','black')
end


% Find change-points with fused Lasso on each profile
disp('Now running flasso on each profile')
figure(2)
clf
for ipro=1:p
    tic
    reslasso = gflassoK(y.profile(:,ipro) , k , o);
    toc
    subplot(nplot,1,ipro);
    plot(y.profile(:,ipro),'o')
    hold on
    ylim=get(gca,'ylim');
    jumps = reslasso.jump';
    line([jumps;jumps],ylim,'linewidth',2,'color','black')
end

disp('Look at the plots')
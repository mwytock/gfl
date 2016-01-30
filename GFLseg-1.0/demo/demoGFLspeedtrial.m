% A script to assess the speed of different variants of group fused Lasso
% with varying signal length, dimension and number of change-points
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert

disp('Now running demoGFLspeedtrial.m')
disp('A speed trial for the different group fused Lasso versions')


% Parameters we vary and their default values

% Length of signals n
nlist = 2.^[4:23];
nn = length(nlist);
ndefault = 1000;

% Number of signals p
plist = 2.^[0:15];
np = length(plist);
pdefault = 10;

% Number of change-points k
klist = 2.^[0:7];
nk = length(klist);
kdefault = 10;


% Number of trials to average over
nTrials = 100;
% Noise level
noiselevel = 0.2;


% Variables to store speed trial results
speedN = zeros(2,nn);
speedP = zeros(2,np);
speedK = zeros(2,nk);

% Initialize the random number generator
s = RandStream.create('mt19937ar','seed',5490);
RandStream.setDefaultStream(s);
    

for trial = 1:nTrials
    disp([num2str(trial),' out of ',num2str(nTrials),' trials'])
    
    % Fix p and k, vary n
    fprintf('N');
    p=pdefault;
    k=kdefault;
    
    counter = 0;
        
        
    % Generate a long random multi-dimensional profile with change points
    nmax = max(nlist);
    y=randomprofile(nmax,k,noiselevel,p);
    
    for n = nlist;
        fprintf('.');
        counter = counter + 1;       

        % Set default weights
        o.weights = defaultweights(n);
        
        % Downsample the profile to length n
        ysmall = y.profile(1:nmax/n:nmax,:);
        
        % Find k change-points with the group fused LARS
        tic;
        reslars = gflars(ysmall,k,o);
        speedN(1,counter) = speedN(1,counter) + toc;
        
        % Find change-points with group fused Lasso at the same lambda
        % value
        tic;
        reslasso = gflasso(ysmall,reslars.lambda(end)-1e-3,o);
        speedN(2,counter) = speedN(2,counter) + toc;
    end
    
    
    % Fix n and k, vary p
    fprintf('P');
    n=ndefault;
    k=kdefault;
    
    % Set default weights
    o.weights = defaultweights(n);
    
    
    counter = 0;

    % Generate a random multi-dimensional profile with change points
    y=randomprofile(n,k,noiselevel,max(plist));

    for p = plist;
        fprintf('.');
        counter = counter + 1;
        
        % Find k change-points with the group fused LARS
        tic;
        reslars = gflars(y.profile(:,1:p),k,o);
        speedP(1,counter) = speedP(1,counter) + toc;
        
        % Find change-points with group fused Lasso at the same lambda
        % value
        tic;
        reslasso = gflasso(y.profile(:,1:p),reslars.lambda(end)-1e-3,o);
        speedP(2,counter) = speedP(2,counter) + toc;
    end
    
    
    
    % Fix n and p, vary k
    fprintf('K');
    n=ndefault;
    p=pdefault;
    
    % Set default weights
    o.weights = defaultweights(n);
    
    % Generate a random multi-dimensional profile with change points
    y=randomprofile(n,max(klist),noiselevel,p);
    
    
    counter = 0;
    
    for k = klist;
        fprintf('.');
        counter = counter + 1;
        
        % Find k change-points with the group fused LARS
        tic;
        reslars = gflars(y.profile,k,o);
        speedK(1,counter) = speedK(1,counter) + toc;
        
        % Find change-points with group fused Lasso at the same lambda values
        tic;
        reslasso = gflasso(y.profile,reslars.lambda(end)-1e-3,o);
        speedK(2,counter) = speedK(2,counter) + toc;
    end
    
    fprintf('\n')
end

% Take average speeds
speedN = speedN/nTrials;
speedP = speedP/nTrials;
speedK = speedK/nTrials;

% Plot performance
figure(1)
clf
loglog(nlist,speedN(1,:))
hold on
loglog(nlist,speedN(2,:),'red')
grid
legend('GFLars','GFLasso','Location','SouthEast')
ylabel('time (s)')
xlabel('n')
hold off


figure(2)
clf
loglog(plist,speedP(1,:))
hold on
loglog(plist,speedP(2,:),'red')
grid
legend('GFLars','GFLasso','Location','SouthEast')
ylabel('time (s)')
xlabel('p')
hold off


figure(3)
clf
loglog(klist,speedK(1,:))
hold on
loglog(klist,speedK(2,:),'red')
grid
legend('GFLars','GFLasso','Location','SouthEast')
ylabel('time (s)')
xlabel('k')
hold off

% % Plots LARS performance
% figure(1);
% subplot(1,3,1)
% plot(nlist,speedN(1,:));
% ylabel('time (s)')
% xlabel('n')
% subplot(1,3,2)
% plot(plist,speedP(1,:));
% ylabel('time (s)')
% xlabel('p')
% subplot(1,3,3)
% plot(klist,speedK(1,:));
% ylabel('time (s)')
% xlabel('k')
% 
% 
% % Plots Lasso performance
% figure(2);
% subplot(1,3,1)
% plot(nlist,speedN(2,:));
% ylabel('time (s)')
% xlabel('n')
% subplot(1,3,2)
% plot(plist,speedP(2,:));
% ylabel('time (s)')
% xlabel('p')
% subplot(1,3,3)
% plot(klist,speedK(2,:));
% ylabel('time (s)')
% xlabel('k')

disp('Look at the plots!')


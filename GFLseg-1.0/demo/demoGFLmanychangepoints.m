%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script estimates the probability of finding 9 jumps correctly
% as a function of the dimension of the signal p and the position of the
% jump u, for fixed noise level. The weighted and unweighted Lars and Lasso are
% compared.
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Now running demoGFLmanychangepoints.m')
disp('A script to estimate the probability of correctly finding several change-point')
disp('If this is too slow, change the "repeats" parameter.')

repeats=100;    % number of repeats to estimate probabilities


%%%%%%%%%%%%%%%%%%
% Data preparation
%%%%%%%%%%%%%%%%%%

% Initialize the random number generator
s = RandStream.create('mt19937ar','seed',5489);
RandStream.setDefaultStream(s);



n=100;   % profile length
beta2=1; % jump variance (jumps are zero-mean Gaussian with variance beta2)
sigma2 = [0.05,0.2,1]';   % noise level
nsigma2 = length(sigma2);

bplist=[10:10:90]'; % Positions of the change-points
nbp = length(bplist);
plist = [1:10:500]; % number of profiles tested
pmax = max(plist);
np = length(plist);

clear acc;
acc.ulars = zeros(nsigma2,np);
acc.wlars = zeros(nsigma2,np);
acc.ulasso = zeros(nsigma2,np);
acc.wlasso = zeros(nsigma2,np);

% Main loop
for r=1:repeats
    fprintf([num2str(r),' out of ',num2str(repeats)]);
    
    % Generate profiles
    Ynonoise = reconstructsignal(n,bplist,randn(nbp,pmax)*sqrt(beta2),ones(n-1,1));
    
    % Loop on the signal level
    for isigma2 = 1:nsigma2
        Y = Ynonoise + randn(n,pmax)*sqrt(sigma2(isigma2));
        

        % Loop on the dimension of the signal
        for ip = 1:np
            p=plist(ip);
            Yp = Y(:,1:p);
            
            % Test unweighted Lars
            o.weights=ones(n-1,nbp);
            r=gflars(Yp,nbp,o); % detect the first breakpoint
            if (sort(r.jump) == bplist)   % and check if is is correct
                acc.ulars(isigma2,ip) = acc.ulars(isigma2,ip) + 1/repeats;
            end
            
            
            % Test weighted Lars
            o.weights = defaultweights(n);
            r=gflars(Yp,nbp,o); % detect the first breakpoint
            if (sort(r.jump) == bplist)   % and check if is is correct
                acc.wlars(isigma2,ip) = acc.wlars(isigma2,ip) + 1/repeats;
            end
            
            % Test unweighted Lasso
            o.weights=ones(n-1,nbp);
            r=gflassoK(Yp,nbp,o); % detect the first breakpoint
            if (sort(r.jump) == bplist)   % and check if is is correct
                acc.ulasso(isigma2,ip) = acc.ulasso(isigma2,ip) + 1/repeats;
            end
            
            % Test weighted Lasso
            o.weights=defaultweights(n);;
            r=gflassoK(Yp,nbp,o); % detect the first breakpoint
            if (sort(r.jump) == bplist)   % and check if is is correct
                acc.wlasso(isigma2,ip) = acc.wlasso(isigma2,ip) + 1/repeats;
            end
            
                    fprintf('.');

        end
    end
    fprintf('\n');
end


%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%
figure
clf
nplot = min(nsigma2,3);
for iplot=1:nplot
    subplot(1,nplot,iplot);
    x1=plot(plist,acc.ulars(iplot,:),'b');
    hold on
    x2=plot(plist,acc.wlars(iplot,:),'r');
    x3=plot(plist,acc.ulasso(iplot,:),'g');
    x4=plot(plist,acc.wlasso(iplot,:),'k');
    hold off
    grid
    axis([0 pmax 0 1])
    xlabel('p')
    ylabel('Accuracy')
    if iplot==1
        legend('U-LARS','W-LARS','U-Lasso','W-Lasso','Location','SouthEast')
    end
end

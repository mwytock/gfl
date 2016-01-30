%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script estimates the probability of finding a single jump correctly
% as a function of the dimension of the signal p and the position of the
% jump u, for fixed noise level. The weighted and unweighted Lars are
% compared (note that they are equivalent to lasso since we only check the
% first jump detected.
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Now running demoGFLsinglechangepoint.m')
disp('A script to estimate the probability of correctly finding a single change-point')
disp('If this is too slow, change the "repeats" parameter.')

repeats=1000;    % number of repeats to estimate probabilities


%%%%%%%%%%%%%%%%%%
% Data preparation
%%%%%%%%%%%%%%%%%%


% Initialize the random number generator
s = RandStream.create('mt19937ar','seed',5489);
RandStream.setDefaultStream(s);


n=100;   % profile length
beta2=1; % common jump size

% Set noise level sigma2 critical for alpha=0.8 in the unweighted case,
% i.e., the correct breakpoint should not be consistently recovered when it
% is at more than 80% of the profile.

alpha=0.8;
sigma2 = n*beta2 * (1-alpha)^2 .* (alpha-1/(2*n)) ./ ( (alpha-0.5-1/(2*n)));


bplist=[50:10:90]; % Where we put the true breakpoints
nbp = length(bplist);
plist = [1:10:500]; % number of profiles tested
pmax = max(plist);
np = length(plist);


acc.unweighted = zeros(nbp,np);
acc.weighted = zeros(nbp,np);
acc.weightedvary = zeros(nbp,np);

%%%%%%%%%%%
% Main loop
%%%%%%%%%%%
for r=1:repeats
    fprintf([num2str(r),' out of ',num2str(repeats)]);
    W = randn(n,pmax)*sqrt(sigma2); % the noise
    for ibp =1:nbp
        bp = bplist(ibp);
        Y=W;
        Y(bp+1:n,:) = Y(bp+1:n,:) + sqrt(beta2); % the profiles
        for ip = 1:np
            p=plist(ip);
            Yp = Y(:,1:p);
            
            % Test unweighted Lars
            o.weights=ones(n-1,1);
            r=gflars(Yp,1,o); % detect the first breakpoint
            if (r.jump(1) == bp)   % and check if it is correct
                acc.unweighted(ibp,ip) = acc.unweighted(ibp,ip) + 1/repeats;
            end
            
            
            % Test weighted Lars
            o.weights = defaultweights(n);
            r=gflars(Yp,1,o); % detect the first breakpoint
            if (r.jump(1) == bp)   % and check if it is correct
                acc.weighted(ibp,ip) = acc.weighted(ibp,ip) + 1/repeats;
            end
            
            
        end
        fprintf('.');
    end
    % Test weighted Lars with varying breakpoint location up to +-2
    for ibp = 1:nbp
        bp = bplist(ibp);
        shift = ones(1,pmax)*bp + round((rand(1,pmax)-1/2)*5);
        Y = W;
        for bb = 1:pmax
            Y(shift(bb)+1:n,bb) = Y(shift(bb)+1:n,bb) + sqrt(beta2); % the profiles
        end
        
        for ip = 1:np
            p=plist(ip);
            Yp = Y(:,1:p);
            
            o.weights = defaultweights(n);
            r=gflars(Yp,1,o); % detect the first breakpoint
            if (abs(r.jump(1)-bp) <= 2)   % and check if it is correct
                acc.weightedvary(ibp,ip) = acc.weightedvary(ibp,ip) + 1/repeats;
            end  
        end
        fprintf('.');

    end
    
    fprintf('\n');
end

%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%
figure
clf
subplot(1,3,1);
plot(plist,acc.unweighted');
grid
axis([0 pmax 0 1])
xlabel('p')
ylabel('Accuracy: unweighted')
legend(strcat('u=',num2str(bplist')),'Location','East')
subplot(1,3,2);
plot(plist,acc.weighted');
grid
axis([0 pmax 0 1])
xlabel('p')
ylabel('Accuracy: weighted')
legend(strcat('u=',num2str(bplist')),'Location','East')
subplot(1,3,3);
plot(plist,acc.weightedvary');
grid
axis([0 pmax 0 1])
xlabel('p')
ylabel('Accuracy: weighted+vary')
legend(strcat('u=',num2str(bplist'),'±2'),'Location','East')

disp('Look at the plot!')
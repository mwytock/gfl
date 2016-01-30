function demoTVL2

%
% Demo TV/L2 solve
%

clear; close all;
path(path,genpath(pwd));


%% generate data -- see nested function below
[I,H,Bn,mu] = genData;
snrBn = snr(Bn,I);
%
%% Run FTVd_v4.0
t = cputime;
out = FTVd_v4(H,Bn,mu,'L2');
t = cputime - t;

%% plot result
figure(1);
subplot(121); imshow(Bn,[]);
title(sprintf('SNR %4.2fdB',snrBn),'fontsize',13);
subplot(122); imshow(out.sol,[]);
title(sprintf('SNR %4.2fdB, CPU %4.2fs, It: %d',snr(out.sol),t,out.itr),'fontsize',13);

figure(2);
subplot(121); plot(1:length(out.snr),out.snr,'b--',1:length(out.snr),out.snr,'r.');
title('SNR history','fontsize',13);
subplot(122); semilogy(1:length(out.f)-1,out.f(2:end),'b--',1:length(out.f)-1,out.f(2:end),'r.');
title('Function values','fontsize',13);
axis([0,inf,min(out.f(2:end)),max(out.f(2:end))]);

fprintf('SNR(Bn) %4.2fdB, SNR(Recovered) %4.2fdB,',snrBn,snr(out.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t,out.itr);

%% nested function
    function [I,H,Bn,mu] = genData
        
        I = double(imread('cameraman.tif'))/255;
        
        [m,n] = size(I);
        H = fspecial('average',15);
        
        sigma = 1.e-3;
        Bn = imfilter(I,H,'circular','conv') + sigma*randn(m,n);
        
        mu = 5.e4;
        
    end

end


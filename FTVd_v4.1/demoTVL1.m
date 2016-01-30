function demoTVL1
%
% Demo TV/L1 solve
%

clear; close all; 
path(path,genpath(pwd));

%% generate data -- see nested function below
[I,H,Bn,mu] = genData;

%% Run FTVd_v4.0
t = cputime;
out = FTVd_v4(H,Bn,mu,'L1');
t = cputime - t;

%% Plot result
figure(1);
subplot(121); imshow(Bn,[]);
title(sprintf('Corruption: %4.0f%%',level*100),'fontsize',13); 
subplot(122); imshow(out.sol,[]);
title(sprintf('SNR %4.2fdB, CPU %4.2fs, It %d',snr(out.sol),t,out.itr),'fontsize',13); 

figure(2);
subplot(121); plot(1:length(out.snr),out.snr,'b--',1:length(out.snr),out.snr,'r.');
title('SNR history','fontsize',13);
subplot(122); semilogy(1:length(out.f)-1,out.f(2:end),'b--',1:length(out.f)-1,out.f(2:end),'r.');
title('Function values','fontsize',13);
axis([0,inf,min(out.f(2:end)),max(out.f(2:end))]);

fprintf('Corruption %2.0d%%, SNR(Recovered) %4.2fdB,',level*100,snr(out.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t,out.itr);

%% nested function 
function [I,H,Bn,mu] = genData
I = double(imread('cameraman.tif'))/255;

% % Trucated kernel (its size is much smaller than that of the true image)
% H = fspecial('average',9);

% % Un-trucated kernel (has the same size with the true image)
[m,n] = size(I);
gauss_std = 3;
wx = exp(-((1:m)-ceil(m/2)).^2/(2*gauss_std^2))';
wy = exp(-((1:n)-ceil(n/2)).^2/(2*gauss_std^2))';
H = wx*wy'; H = H/sum(sum(H));

level = .6;
Bn = imfilter(I,H,'circular','conv');
Bn = imnoise(Bn,'salt & pepper',level);
snr(Bn,I);

% suggested mu (not too bad)
switch level
    case 0.10; mu = 60;
    case 0.20; mu = 40;
    case 0.30; mu = 13;
    case 0.40; mu = 10;
    case 0.50; mu = 8;
    case 0.60; mu = 4;
    otherwise
        mu = input('Input mu:');
end

end

end

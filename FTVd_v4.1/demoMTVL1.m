function demoMTVL1
%
% MTV/L1 solve
%

clear; close all;
path(path,genpath(pwd));

%% generate data -- see nested function below
[I,H,Bn,mu] = genData;

%% run FTVd_v4.0
t = cputime;
out = FTVd_v4(H,Bn,mu,'L1');
t = cputime - t;

%% plot result
subplot(121);
imshow(Bn,[]);
title(sprintf('Corruption: %4.0f%%',level*100),'fontsize',13); 

subplot(122);
imshow(out.sol,[]);
title(sprintf('SNR %4.2fdB, CPU %4.2fs, It %d, ', ...
    snr(out.sol),t,out.itr),'fontsize',13);

fprintf('Corruption %2.0d%%, SNR(Recovered) %4.2fdB,',level*100,snr(out.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t,out.itr);

%% nested function
    function [I,H,Bn,mu] = genData
        
        % original image
        I = imread('lena256.png');
        I = double(I)/255;
        
        % parameters for 9 kernels
        len = 21; theta = 45;
        len2 = 41; theta2 = 90;
        len3 = 61; theta3 = 135;
        
        hsize = 11; sigma = 9;
        hsize2 = 21; sigma2 = 11;
        hsize3 = 31; sigma3 = 13;
        
        ahsize = 13;
        ahsize2 = 15;
        ahsize3 = 17;
        
        % 9 kernels
        Htmp = cell(1,9);
        Htmp{1} = fspecial('motion',len,theta);
        Htmp{2} = fspecial('motion',len2,theta2);
        Htmp{3} = fspecial('motion',len3,theta3);
        Htmp{4} = fspecial('gaussian',hsize,sigma);
        Htmp{5} = fspecial('gaussian',hsize2,sigma2);
        Htmp{6} = fspecial('gaussian',hsize3,sigma3);
        Htmp{7} = fspecial('average',[ahsize, ahsize]);
        Htmp{8} = fspecial('average',[ahsize2,ahsize2]);
        Htmp{9} = fspecial('average',[ahsize3,ahsize3]);
        
        % weights assigned to each kernel
        r = [ .7 .15 .15;
            .1 .8 .1;
            .2 .2 .6];
        
        % kernel
        H = cell(3,3);
        R = randperm(9);
        H{1,1} = Htmp{R(1)}; H{1,2} = Htmp{R(2)}; H{1,3} = Htmp{R(3)};
        H{2,1} = Htmp{R(4)}; H{2,2} = Htmp{R(5)}; H{2,3} = Htmp{R(6)};
        H{3,1} = Htmp{R(7)}; H{3,2} = Htmp{R(8)}; H{3,3} = Htmp{R(9)};
        
        H{1,1} = r(1,1).*H{1,1}; H{1,2} = r(1,2).*H{1,2}; H{1,3} = r(1,3).*H{1,3};
        H{2,1} = r(2,1).*H{2,1}; H{2,2} = r(2,2).*H{2,2}; H{2,3} = r(2,3).*H{2,3};
        H{3,1} = r(3,1).*H{3,1}; H{3,2} = r(3,2).*H{3,2}; H{3,3} = r(3,3).*H{3,3};
        
        % Observation
        level = .6;
        Bn = imfilter33(I,H,'conv');
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

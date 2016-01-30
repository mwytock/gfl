% Run timing tests on real and synthetic data

randn('seed',1);
rand('seed',1);

% Compare to FTVd on a single image
img = 'lena.bmp';   w = 256; k_max = w;
img = 'colors.png'; w = 256; k_max = 100;
I = double(imresize(imread(img), [w nan]))/256;
Y0 = problemDenoise(img, w);

F = permute(Y0, [2 3 1]);
snr0 = snr(F,I)
imwrite(F, ['noise_' img]);

% Run FTVd
% Identity convolution
% H = {1, 0, 0
%      0, 1, 0
%      0, 0, 1};
% mu = 10;
% opts = struct;
% opts.print = 1;
% out1 = ADM2MTVL2(H,F,mu,opts);

% Run Dykstra + ASPN
lambda = 0.2;
lam = lambda*ones(w-1,1);
Y = Y0;
opts = options('aspn');
opts.verbose = 0;
opts.tol = 1e-6;
opts.k_max = k_max;
out2 = denoise3d(Y,lam,lam,opts);

% Run ADMM with lattice operator for D
D = lattice_operator(w,w);
Y = Y0(:,:)';
opts = options('adm');
opts.reltol = 1e-6;
opts.abstol = 1e-6;
opts.rho = 2;
opts.maxit = 1000;
out3 = admmDSolver(Y,D,lambda,opts);
save(['denoise_' img '.mat'], 'w', 'out2', 'out3');
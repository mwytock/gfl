function opt=options(alg)
% OPTIONS -- generate options struct for specified algo

   if ~exist('alg', 'var')
      alg='pn';
   end

   opt.maxit=1e4;                     % number of iterations
   opt.stats=1;                         % compute stats like objval, dgap
   opt.alg=alg;
   opt.max_time=300;
   opt.verbose=1;

   switch alg
     case 'pg'
       opt.solver='Projected Gradient';
     case 'og'
       opt.solver='Accelerated Projected Gradient';
     case 'spg'
       opt.solver='Spectral Projected Gradient';
     case 'dr'
       opt.solver='Douglas-Rachford';
     case 'lb'
       opt.solver='LFBGS-B';
       opt.factr=1e4;
       opt.pgtol=1e-8;
       opt.m=10;
       opt.maxIts=opt.maxit;
     case 'aspn'
       opt.solver='Active Set Projected Newton';
       opt.verbose = 1;
       opt.max_iters = opt.maxit;
       opt.max_ls_iters = 1000;
       opt.eps = 1e-8;
       opt.tol = 1e-6;
       opt.beta = 0.5;
       opt.k_max = 500;
       opt.ls_quadratic = 1;
     case 'pn'
       opt.solver='Projected Newton';
       opt.verbose = 1;
       opt.max_iters = opt.maxit;
       opt.max_ls_iters = 1000;
       opt.eps = 1e-8;
       opt.tol = 1e-4;
       opt.beta = 0.5;
       opt.k_max = 500;
       opt.ls_quadratic = 1;
     case 'sbb'
       opt.solver='Subspace BB';
       opt.tolg=1e-8;
       opt.truex=0;
       opt.verbose=1;
     case 'plb'
       opt.solver='Projected LBFGS-B';
       opt.tolinfg=1e-8;       % convergence check: norm(gradient, 'inf') <= opt.tolinfg
       opt.usetolginf = 1;     % use as default convg. criterion
       opt.tolo       = 1e-12; % stopping criterion based on change in obj
       opt.maxmem = 11;           % maximum number of LBFGS memory vectors
       opt.max_func_evals = 30;
       opt.beta = 0.0498;
       opt.sigma = 0.298;
       opt.eps=1e-8;
       opt.grad=[];
     case 'gflasso'
       opt.solver='GFLseg';
       opt.tol = 1e-8;
       opt.verbose = 1;
     case 'trip'
       opt.maxit=20;
       opt.tolo       = 1e-9;   % stopping criterion based on change in obj
       opt.tolx       = 1e-6;   % stopping criterion based on change in x
       opt.tolro      = 0;      % stopping criterion based on relative change in obj
       opt.maxt       = inf;    % stopping criterion based on running time
       % semi-advanced params
       opt.maxnull    = 50;    % maximum number of null iteration
       opt.alpmax     = 1e+8;  % maximum BB step-size
       opt.alpmin     = 1e-8;  % minimum BB step-size
       opt.tollcp     = 1e-2;   % stopping criterion based on LCP

       % advanced, better not change unless you know what you are doing!
       opt.del0       = 1000;    % initial trust-region
       opt.delmax     = 1e+5;   % maximum trust-region
       opt.delmin     = 1e-5;   % minimum trust-region
       opt.toldel     = 1e-9;   % stopping criterion based on the radius of TR
       opt.eta1       = 0.01;   % constant to determine model accuracy
       opt.eta2       = 0.9;    % constant to determine model accuracy
       opt.gam1       = 0.5;    % -
       opt.gam2       = 0.5;    % -
       opt.gam3       = 2;      % - TR update constants.
       opt.beta       = 5e-1;   % weight to determin sufficient model reduction
       opt.tau        = .37;    % step-reduction rate in the line-search
       opt.std        = 0;      % standardize the data
     case 'adm'
       opt.solver='ADMM';
       opt.rho=1;
       opt.abstol=1e-4;
       opt.reltol=1e-2;
     otherwise
       opt=[];
       fprintf('ERROR: options(%s) undefined\n',alg);
   end
end
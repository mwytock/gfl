function out = ADM2MTVL2(H,F,mu,opts)
% out = ALD2MTVL2(H,F,mu,opts);
%
% Alternating Directions Method (ADM) applied to MTV/L2.
%
% Suppose the data accuquisition model is given by: F = K*Xbar + Noise,
% where Xbar is an original color image, K is a cross-channel convolution
% matrix, Noise is additive, and F is a blurry and noisy observation. To
% recover Xbar from F and K, we solve MTV/L2 model
%
% ***     min_X \sum_i ||kron(I3,Di)*X|| + mu/2*||K*X - F||^2      ***
%
% Inputs:
%         H  ---  convolution kernel representing K
%         F  ---  blurry and noisy observation
%         mu ---  model prameter (must be provided by user)
%         opts --- a structure containing algorithm parameters {default}
%                 * opst.beta    : a positive constant {10}
%                 * opst.gamma   : a constant in (0,1.618] {1.618}
%                 * opst.maxitr  : maximum iteration number {500}
%                 * opst.relchg  : a small positive parameter which controls
%                                  stopping rule of the code. When the
%                                  relative change of X is less than
%                                  opts.relchg, then the code stops. {1.e-3}
%                 * opts.print   : print inter results or not {1}
%
% Outputs:
%         out --- a structure contains the following fields
%                * out.snr   : SNR values at each iteration
%                * out.relchg: the history of relative change in X
%                * out.sol   : numerical solution obtained by this code
%                * out.itr   : number of iterations used
%

% Copyright (c), May, 2009
%       Junfeng Yang, Dept. Math., Nanjing Univiversity
%       Wotao Yin,    Dept. CAAM, Rice University
%       Yin Zhang,    Dept. CAAM, Rice University


[m n d3] = size(F);

if d3 == 1
    error('Call ADM2TVL2 please!');
end

% default parameter settings
if nargin < 4; opts = []; end
opts = getopts(opts);

% initialization
X = F;
Lam1 = zeros(m,n,d3);
Lam2 = Lam1;

% algorithm parameters
beta = opts.beta;
gamma = opts.gamma;
print = opts.print;

% compute constant quantities
C = getC;

% define finite difference operator for color image
[D,Dt] = defMDDt;

% finite diff.
[D1X,D2X] = D(X);

% keep record
out.snr = [];
out.relchg = [];
out.time0 = cputime;
out.time = [];

%% Main loop
for ii = 1:opts.maxitr

    % ==================
    %   Shrinkage Step
    % ==================
    Z1 = D1X + Lam1/beta;
    Z2 = D2X + Lam2/beta;

    V3 = sum(Z1.^2 + Z2.^2,d3);
    V3 = sqrt(V3);
    V = zeros(m,n,d3);
    for jj = 1:d3;
        V(:,:,jj) = V3;
    end
    clear V3
    V(V==0) = 1;
    V = max(V - 1/beta, 0)./V;
    Y1 = Z1.*V;
    Y2 = Z2.*V;

    % ==================
    %     X-subprolem
    % ==================
    Xp = X;
    X = (mu*C.KtF - Dt(Lam1,Lam2))/beta + Dt(Y1,Y2);
    X = getFX;
    X = real(ifft2(X));

    %% keep record
    relchg = norm(Xp(:) - X(:))/norm(Xp(:));
    out.relchg = [out.relchg; relchg];
    snrX = snr(X);
    out.snr = [out.snr; snrX];
    out.time = [out.time; cputime];

    if print
        fprintf('Iter: %d, snrX: %4.2f, relchg: %4.2e\n',ii,snrX,relchg);
    end

    %% check stopping rule
    if relchg < opts.relchg
        out.sol = X;
        out.itr = ii;
        return
    end

    % ==================
    %    Update Lam
    % ==================
    [D1X,D2X] = D(X);
    Lam1 = Lam1 - gamma*beta*(Y1 - D1X);
    Lam2 = Lam2 - gamma*beta*(Y2 - D2X);

end
out.sol = X;
out.itr = ii;
out.exit = 'Exist Normally';
if ii == opts.mxitr
    out.exit = 'Maximum iteration reached!';
end

    function opts = getopts(opts)

        if ~isfield(opts,'maxitr')
            opts.maxitr = 500;
        end
        if ~isfield(opts,'beta')
            opts.beta = 10;
        end
        if ~isfield(opts,'gamma')
            opts.gamma = 1.618;
        end
        if ~isfield(opts,'relchg')
            opts.relchg = 1.e-3;
        end
        if ~isfield(opts,'print')
            opts.print = 0;
        end
    end

%% nested functions
    function C = getC

        % eigenvalues of DtD
        eigsDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;

        % Index
        J1 = 1:m; J2 = (m+1):(2*m); J3 = (2*m+1):(3*m);
        K1 = 1:n; K2 = (n+1):(2*n); K3 = (2*n+1):(3*n);

        % eigenvalues of each K_ij {i,j = 1,2,3}
        eigsK = zeros(3*m,3*n);
        eigsK(J1,K1) = psf2otf(H{1,1},[m,n]);
        eigsK(J1,K2) = psf2otf(H{1,2},[m,n]);
        eigsK(J1,K3) = psf2otf(H{1,3},[m,n]);
        eigsK(J2,K1) = psf2otf(H{2,1},[m,n]);
        eigsK(J2,K2) = psf2otf(H{2,2},[m,n]);
        eigsK(J2,K3) = psf2otf(H{2,3},[m,n]);
        eigsK(J3,K1) = psf2otf(H{3,1},[m,n]);
        eigsK(J3,K2) = psf2otf(H{3,2},[m,n]);
        eigsK(J3,K3) = psf2otf(H{3,3},[m,n]);

        % eigenvalues of each block in KtK
        eigsKtK = zeros(3*m,3*n);
        eigsKtK(J1,K1) = abs(eigsK(J1,K1)).^2 + abs(eigsK(J2,K1)).^2 + abs(eigsK(J3,K1)).^2;
        eigsKtK(J1,K2) = conj(eigsK(J1,K1)).*eigsK(J1,K2) + conj(eigsK(J2,K1)).*eigsK(J2,K2) + conj(eigsK(J3,K1)).*eigsK(J3,K2);
        eigsKtK(J1,K3) = conj(eigsK(J1,K1)).*eigsK(J1,K3) + conj(eigsK(J2,K1)).*eigsK(J2,K3) + conj(eigsK(J3,K1)).*eigsK(J3,K3);

        eigsKtK(J2,K1) = conj(eigsK(J1,K2)).*eigsK(J1,K1) + conj(eigsK(J2,K2)).*eigsK(J2,K1) + conj(eigsK(J3,K2)).*eigsK(J3,K1);
        eigsKtK(J2,K2) = abs(eigsK(J1,K2)).^2 + abs(eigsK(J2,K2)).^2 + abs(eigsK(J3,K2)).^2;
        eigsKtK(J2,K3) = conj(eigsK(J1,K2)).*eigsK(J1,K3) + conj(eigsK(J2,K2)).*eigsK(J2,K3) + conj(eigsK(J3,K2)).*eigsK(J3,K3);

        eigsKtK(J3,K1) = conj(eigsK(J1,K3)).*eigsK(J1,K1) + conj(eigsK(J2,K3)).*eigsK(J2,K1) + conj(eigsK(J3,K3)).*eigsK(J3,K1);
        eigsKtK(J3,K2) = conj(eigsK(J1,K3)).*eigsK(J1,K2) + conj(eigsK(J2,K3)).*eigsK(J2,K2) + conj(eigsK(J3,K3)).*eigsK(J3,K2);
        eigsKtK(J3,K3) = abs(eigsK(J1,K3)).^2 + abs(eigsK(J2,K3)).^2 + abs(eigsK(J3,K3)).^2;

        % eigenvalues of each block in DtD + (mu/beta)*KtK
        C.eigs = (mu / beta) * eigsKtK;
        C.eigs(J1,K1) = C.eigs(J1,K1) + eigsDtD;
        C.eigs(J2,K2) = C.eigs(J2,K2) + eigsDtD;
        C.eigs(J3,K3) = C.eigs(J3,K3) + eigsDtD;

        % clear variables to save space
        clear eigsDtD eigsKtK eigsK

        % KtF
        C.KtF = imfilter33(F,H,'corr');

        % The following lines are used to diagonalize the block
        % diagonalized matrix DtD + (mu/beta) KtK, which is useful in
        % computing FX by block Gaussian elimination at each iteration.
        C.Tmp = cell(1,6);
        Tmp1 = C.eigs(J2,K1)./C.eigs(J1,K1);
        Tmp2 = C.eigs(J3,K1)./C.eigs(J1,K1);
        C.eigs(J2,K2) = C.eigs(J2,K2) - C.eigs(J1,K2).*Tmp1;
        C.eigs(J3,K2) = C.eigs(J3,K2) - C.eigs(J1,K2).*Tmp2;
        C.eigs(J2,K3) = C.eigs(J2,K3) - C.eigs(J1,K3).*Tmp1;
        C.eigs(J3,K3) = C.eigs(J3,K3) - C.eigs(J1,K3).*Tmp2;
        Tmp3 = C.eigs(J3,K2)./C.eigs(J2,K2);
        Tmp4 = C.eigs(J1,K2)./C.eigs(J2,K2);
        C.eigs(J3,K3) = C.eigs(J3,K3) - C.eigs(J2,K3).*Tmp3;
        C.eigs(J1,K3) = C.eigs(J1,K3) - C.eigs(J2,K3).*Tmp4;
        Tmp5 = C.eigs(J2,K3)./C.eigs(J3,K3);
        Tmp6 = C.eigs(J1,K3)./C.eigs(J3,K3);

        C.Tmp{1,1} = Tmp1;
        C.Tmp{1,2} = Tmp2;
        C.Tmp{1,3} = Tmp3;
        C.Tmp{1,4} = Tmp4;
        C.Tmp{1,5} = Tmp5;
        C.Tmp{1,6} = Tmp6;
    end

    function FX = getFX
        % Compute Fourier transform of the new iterate
        FX = fft2(X);
        FX(:,:,2) = FX(:,:,2) - C.Tmp{1,1}.*FX(:,:,1);
        FX(:,:,3) = FX(:,:,3) - C.Tmp{1,2}.*FX(:,:,1);
        FX(:,:,3) = FX(:,:,3) - C.Tmp{1,3}.*FX(:,:,2);
        FX(:,:,1) = FX(:,:,1) - C.Tmp{1,4}.*FX(:,:,2);
        FX(:,:,2) = FX(:,:,2) - C.Tmp{1,5}.*FX(:,:,3);
        FX(:,:,1) = FX(:,:,1) - C.Tmp{1,6}.*FX(:,:,3);

        J1 = 1:m; J2 = (m+1):(2*m); J3 = (2*m+1):(3*m);
        K1 = 1:n; K2 = (n+1):(2*n); K3 = (2*n+1):(3*n);
        FX(:,:,1) = FX(:,:,1)./C.eigs(J1,K1);
        FX(:,:,2) = FX(:,:,2)./C.eigs(J2,K2);
        FX(:,:,3) = FX(:,:,3)./C.eigs(J3,K3);
    end

    function [D,Dt] = defMDDt
        % finite difference operator D
        % and its transpose operator
        % for Multichannel images

        D = @(U) ForwardD(U);
        Dt = @(X,Y) Dive(X,Y);

        function [Dux,Duy] = ForwardD(U)
            Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
            Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
        end

        function DtXY = Dive(X,Y)
            DtXY = [X(:,end,:) - X(:,1,:), -diff(X,1,2)];
            DtXY = DtXY + [Y(end,:,:) - Y(1,:,:); -diff(Y,1,1)];
        end

    end

end

function out = ADM2TVL1(H,F,mu,opts)
% out = ADM2TVL1(H,F,mu,opts);
%
% Alternating Directions Method (ADM) applied to TV/L1.
% 
% Suppose the data accuquisition model is given by: F = N(K*Xbar),
% where Xbar is an original image, K is a convolution matrix, N(cdot)
% represents a procedure of impulsive noise corruption, and F is 
% a blurry and noisy observation. To recover Xbar from F and K, 
% we solve TV/L1 model
% 
% ***     min_X \sum_i ||Di*X|| + mu*||K*X - F||_1      ***
% 
% Inputs:
%         H  ---  convolution kernel representing K
%         F  ---  blurry and noisy observation
%         mu ---  model prameter (must be provided by user)
%         opts --- a structure containing algorithm parameters {default}
%                 * opst.beta1   : a positive constant {5}
%                 * opst.beta2   : a positive constant {100}
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
%                * out.f     : function valuse at each itertion
%                * out.tv    : TV values at each iteration
%                * out.fid   : fidelity values at each iteration
%                * out.relchg: the history of relative change in X
%                * out.sol   : numerical solution obtained by this code
%                * out.itr   : number of iterations used
% 

% Copyright (c), May, 2009
%       Junfeng Yang, Dept. Math., Nanjing Univiversity
%       Wotao Yin,    Dept. CAAM, Rice University
%       Yin Zhang,    Dept. CAAM, Rice University

[m n d3] = size(F);

if d3 == 3
    error('Call ADM2MTVL1 please!');
end

if nargin < 4; opts = []; end
opts = getopts(opts);

C = getC(F,H); 
[D,Dt] = defDDt;

% initialization
X = F;
Lam1 = zeros(m,n);
Lam2 = Lam1;
Lam3 = Lam1;
beta1 = opts.beta1;
beta2 = opts.beta2;
gamma = opts.gamma;
print = opts.print;

Denom = C.eigsDtD + beta2/beta1 * C.eigsKtK;

% finite diff
[D1X,D2X] = D(X);
KXF = real( ifft2(C.eigsK .* fft2(X)) ) - F;

[tv,fid,f] = fval(D1X,D2X,KXF,mu);

out.snr = [];
out.nrmdX = [];
out.relchg = [];
out.f = f;
out.tv = tv;
out.fid = fid;

%% Main loop
for ii = 1:opts.maxitr

    V1 = D1X + Lam1/beta1;
    V2 = D2X + Lam2/beta1;
    V3 = KXF + Lam3/beta2;

    V = V1.^2 + V2.^2;
    V = sqrt(V);
    V(V==0) = 1;

    % ==================
    %   Shrinkage Step
    % ==================
    V = max(V - 1/beta1, 0)./V;
    Y1 = V1.*V;
    Y2 = V2.*V;
    Z = max(abs(V3) - mu/beta2, 0).*sign(V3);

    % ==================
    %     X-subprolem
    % ==================
    Xp = X;
    Temp = (beta2*Z - Lam3)/beta1; 
    Temp = real( ifft2(C.eigsKt .* fft2(Temp)) );
    
    X = Dt(Y1 - Lam1/beta1,Y2 - Lam2/beta1) + Temp + beta2/beta1*C.KtF;
    X = fft2(X)./Denom;
    X = real(ifft2(X));

    snrX = snr(X);
    out.snr = [out.snr; snrX];
    relchg = norm(X - Xp,'fro')/norm(X,'fro');
    out.relchg = [out.relchg; relchg];

    if print 
        fprintf('Iter: %d, snrX: %4.2f, relchg: %4.2e\n',ii,snrX,relchg);
    end

    % ====================
    % Check stopping rule
    % ====================
    if relchg < opts.relchg
        out.sol = X;
        out.itr = ii;
        [D1X,D2X] = D(X);
        KXF = real( ifft2(C.eigsK .* fft2(X)) ) - F;
        [tv,fid,f] = fval(D1X,D2X,KXF,mu);
        out.f = [out.f; f];
        out.tv = [out.tv; tv];
        out.fid = [out.fid; fid];
        return
    end

    % finite diff.
    [D1X,D2X] = D(X);
    KXF = real( ifft2(C.eigsK .* fft2(X)) ) - F;
    [tv,fid,f] = fval(D1X,D2X,KXF,mu);
    out.f = [out.f; f];
    out.tv = [out.tv; tv];
    out.fid = [out.fid; fid];

    % ==================
    %    Update Lam
    % ==================
    Lam1 = Lam1 - gamma*beta1*(Y1 - D1X);
    Lam2 = Lam2 - gamma*beta1*(Y2 - D2X);
    Lam3 = Lam3 - gamma*beta2*(Z - KXF);

end
out.sol = X;
out.itr = ii;
out.exit = 'Exist Normally';
if ii == opts.mxitr
    out.exit = 'Maximum iteration reached!';
end

    
%% ------------------SUBFUNCTION-----------------------------
function opts = getopts(opts)

if ~isfield(opts,'maxitr')
    opts.maxitr = 500;
end
if ~isfield(opts,'beta1')
    opts.beta1 = 5;
end
if ~isfield(opts,'beta2')
    opts.beta2 = 20;
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

%% ------------------SUBFUNCTION-----------------------------
function C = getC(F,H)

[m,n] = size(F);

C.eigsDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;
C.eigsK = psf2otf(H,[m,n]);
C.eigsKt = conj(C.eigsK);
C.eigsKtK = abs(C.eigsK).^2;

C.KtF = real(ifft2(C.eigsKt .* fft2(F)));

%% ------------------SUBFUNCTION-----------------------------
function [tv,fid,f] = fval(D1X,D2X,KXF,mu)

tv =  sum(sum(sqrt(D1X.^2 + D2X.^2)));
fid =  sum(abs(KXF(:)));
f = tv + mu * fid;

function [D,Dt] = defDDt

D = @(U) ForwardD(U);
Dt = @(X,Y) Dive(X,Y);

function [Dux,Duy] = ForwardD(U)

Dux = [diff(U,1,2), U(:,1) - U(:,end)];
Duy = [diff(U,1,1); U(1,:) - U(end,:)];

function DtXY = Dive(X,Y)

DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];

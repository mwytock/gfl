function out = FTVd_v4(H,Bn,mu,str,opts)
% out = FTVd_v4(H,Bn,mu,str,opts)
%
% FTVd_v4 solves TV models (TV/L2, TV/L1 and their color variants)
% by Alternating Direction Method (ADM).
%

% Copyright (c), May, 2009
%       Junfeng Yang, Dept. Math., Nanjing Univiversity
%       Wotao Yin,    Dept. CAAM, Rice University
%       Yin Zhang,    Dept. CAAM, Rice University

if nargin < 4
    error('Not enough inputs!');
end

if ~strcmp(str,'L2') && ~strcmp(str,'l2') ...
        && ~strcmp(str,'L1') && ~strcmp(str,'l1')
    error('Model is not specified correctly! See "str".');
end

[m,n,d3] = size(Bn);

if d3 ~= 1 && d3 ~= 3
    error('FTVd applies to gray and RGB color images only!');
end


if nargin < 5
    opts = [];
end


disp('FTVd_v4 is running, please wait ...');


if d3 == 1 % grayscale image
    if strcmp(str,'L2') || strcmp(str,'l2')
        out = ADM2TVL2(H,Bn,mu,opts);
    else
        out = ADM2TVL1(H,Bn,mu,opts);
    end
else % RGB image
    if strcmp(str,'L2') || strcmp(str,'l2')
        out = ADM2MTVL2(H,Bn,mu,opts);
    else
        out = ADM2MTVL1(H,Bn,mu,opts);
    end
end

disp('Done!');



function B = imfilter33(U,H,s)
% B = imfilter33(U,H,s);
% Compute Cross-channel convolution H*U
%
% Inputs:
%       -- s:    'corr', 'conv'
%       -- U:    [ Ur; Ug; Ub ];
%       -- H:  a 3 x 3 cell array carrying cross-channel blurring kernels:
%               [ H{1,1} H{1,2} H{1,3}; 
%                H{2,1} H{2,2} H{2,3}; 
%                H{3,1} H{3,2} H{3,3} ];
%
% Output -- B = H*U
%

[m,n,d3] = size(U);

if d3 ~= 3;
    error(['Underlying image does not have 3 channels! ', ...
        'For siglechannel image, please use imfilter.m']);
elseif ~iscell(H)
    Htmp = H;
    H = cell(3,3);
    H{1,1} = Htmp;
    H{2,2} = Htmp;
    H{3,3} = Htmp;
else
    [ah,bh] = size(H);
    if ah ~= 3 || bh ~= 3
        error('H is not a 3 x 3 cell array');
    end
end
    

B = zeros(m,n,d3);
if nargin < 3
    error('Input not enough for imfilter33!');
end


if  strcmp(s,'corr') 
    B(:,:,1) = imfilter(U(:,:,1),H{1,1},'circular','corr') + ...
               imfilter(U(:,:,2),H{2,1},'circular','corr') + ...
               imfilter(U(:,:,3),H{3,1},'circular','corr');

    B(:,:,2) = imfilter(U(:,:,1),H{1,2},'circular','corr') + ...
               imfilter(U(:,:,2),H{2,2},'circular','corr') + ...
               imfilter(U(:,:,3),H{3,2},'circular','corr');

    B(:,:,3) = imfilter(U(:,:,1),H{1,3},'circular','corr') + ...
               imfilter(U(:,:,2),H{2,3},'circular','corr') + ...
               imfilter(U(:,:,3),H{3,3},'circular','corr');
elseif strcmp(s,'conv')
    B(:,:,1) = imfilter(U(:,:,1),H{1,1},'circular','conv') + ...
               imfilter(U(:,:,2),H{1,2},'circular','conv') + ...
               imfilter(U(:,:,3),H{1,3},'circular','conv');

    B(:,:,2) = imfilter(U(:,:,1),H{2,1},'circular','conv') + ...
               imfilter(U(:,:,2),H{2,2},'circular','conv') + ...
               imfilter(U(:,:,3),H{2,3},'circular','conv');

    B(:,:,3) = imfilter(U(:,:,1),H{3,1},'circular','conv') + ...
               imfilter(U(:,:,2),H{3,2},'circular','conv') + ...
           imfilter(U(:,:,3),H{3,3},'circular','conv');
else
    error('convolution of correlation? please check imfilter33.m');
end


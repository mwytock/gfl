function x = snr(sig, ref)
% x = snr(sig, ref)
% snr -- Compute Signal-to-Noise Ratio for images
%
% Usage:
%       x = snr(sig, ref)  -- 1st time or
%       x = snr(sig)       -- afterwards
%
% Input:
%       sig         Modified image
%       ref         Reference image
%  
% Output:
%       x           SNR value

persistent ref_save;

if nargin == 2; ref_save = ref; end;

mse = mean((ref_save(:)-sig(:)).^2);
dv = var(ref_save(:),1);
x = 10*log10(dv/mse);

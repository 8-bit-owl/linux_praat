function x = createBandlimitedNoise(L,fs,noiseBandwidth,noiseFilterOrder)
% This function creates band-limited Gaussian random noise of length L samples.  If noiseBandwith is
% greater than fs/2, no filtering is applied.
% Input:
%    L -- integer scalar, samples, length of noise vector
%    fs -- real scalar, Hz, sampling rate of noise vector
%    noiseBandwidth -- real scalar, Hz, base bandwidth of noise vector
%    noiseFilterOrder -- integer scalar, order of Butterworth filter to use for bandlimiting (Note:
%       second-order system filtering is applied to avoid instabilites of high-order Butterworth
%       filters).
% Output:
%    x -- Lx1 real vector, Gaussian random signal, bandlimited if necessary, drawn from zero mean,
%    unity variance distribution using randn() function.

% Mark Skowronski, October 18, 2013

% Create noise vector:
x = randn(L,1);

% Filter if necessary:
if noiseBandwidth<(fs/2),
   % Design Butterworth lowpass filter:
   [z,p,k] = butter(noiseFilterOrder,noiseBandwidth/(fs/2)); % pole-zero representation
   
   % Filter x using second-order system filtering:
   x = xFilter(x,z,p,k);
end;

return;
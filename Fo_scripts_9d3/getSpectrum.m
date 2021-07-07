function [X,fX,N] = getSpectrum(x,fs,N,windowType)
% This function calculates the power spectral density of x and FFT bin frequencies.
% Input [defaults]:
%    x -- Lx1 real vector, time domain signal
%    fs -- real scalar, Hz, sampling rate of x
%    N -- integer scalar, FFT size of x (if N<length(x), N set to [length(x)])
%    windowType -- string of FFT window name: ['hamming'],'hanning','blackman','rectangular'
% Output:
%    X -- Qx1 real vector, left half of power spectral density of x, Q=N/2+1 (N even) or (N-1)/2+1
%         (N odd)
%    fX -- Qx1 real vector, Hz, bin frequencies of X
%    N -- integer scalar, max(N input,length(x)), actually used in FFT construction
%
% Power spectral density of x = abs(fft(x.*w,N)).^2/(L*U), U = mean-squared value of window w

% Get window function:
L = length(x); % samples, window length
w = getWindowFunction(windowType,L);
U = mean(w.^2);

% Construct FFT of x:
N = max(N,length(x)); % ensure FFT size is at least length(x)
xFFT = fft(x(:).*w,N);

% Truncate FFT to left half of spectrum, construct bin frequencies:
if mod(N,2)==0, % even FFT length
   Q = N/2+1; % includes DC, lower half, and Nyquist frequencies
else % odd FFT length
   Q = (N-1)/2+1; % includes DC and lower half frequencies
end;
xFFT = xFFT(1:Q);
fX = [0:Q-1]'*fs/N; % Hz, frequency bins of X, COLUMN vector

% Convert FFT to power spectral density:
X = (real(xFFT).^2+imag(xFFT).^2)/(L*U); % faster than abs(xFFT).^2

return;

function w = getWindowFunction(windowType,L)
% This function returns a COLUMN vector of a window function useful for spectral analysis.

switch lower(windowType)
   case 'blackman'
      w = blackman(L);
   case 'hanning'
      w = hanning(L);
   case 'rectangular'
      w = rectwin(L);
   otherwise % Hamming window is default
      w = hamming(L);
end;

return;


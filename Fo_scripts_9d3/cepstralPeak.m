function [CP,F0,parameters,c] = cepstralPeak(x,fs,parameters)
% This function calculates the cepstral peak for the input.
% Input:
%    x -- Nx1 real vector, data vector to process
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters [defaults]
%       .windowType -- ['rectangular'],'Hamming','Hanning','Blackman', temporal window applied to x
%       .FFTsize -- real scalar, samples, size of FFT to apply to x [1024]
%       .clipLimit -- real scalar, dB, spectrum floor set at spectral peak minus clipLimit, set to
%          [] to ignore spectral floor clipping [100]
%       .IFFTfactor -- real scalar, IFFT size is FFTsize*IFFTfactor, for cepstral interpolation,
%          should be factor of 2 for IFFT efficiency [8]
%       .F0range -- 1x2 real vector, Hz, range of F0 over which to search for peak [50 400]
% Output:
%    CP -- real scalar, cepstral peak
%    F0 -- real scalar, Hz, frequency at CP
%    parameters -- same as input (or default values if empty)
%    c -- Mx1 real vector, real cepstrum of x, M=FFTsize*IFFTfactor

% Check parameters:
if nargin<3,
   parameters(1).windowType = 'rectangular';
   parameters(1).FFTsize = 1024;
   parameters(1).clipLimit = 100; % dB
   parameters(1).IFFTfactor = 8;
   parameters(1).F0range = [50 400]; % Hz
end;

% Check x:
if min(size(x))>1, % multi-dimensional x
   % Display error:
   error('x must be COLUMN vector.');
else
   x = x(:); % Ensure x is COLUMN vector
end;

% Check if FFTsize > length(x), warn if not:
if size(x,1)>parameters.FFTsize,
   FFTsizeRec = 2^(nextpow2(size(x,1))+1); % recommended FFTsize
   warning(['FFTsize shorter than x length, resulting in x truncation.',...
      ' Increase FFTsize to eliminate truncation, recommend: FFTsize = ',num2str(FFTsizeRec),'.']);
end;

% Create window:
w = createWindow(parameters.windowType,length(x));

% Apply window, make log magnitude spectrum:
X = fft(x.*w,parameters.FFTsize);
Xlog = log(abs(X));

% Apply spectral floor if necessary:
if ~isempty(parameters.clipLimit),
   Xfloor = max(Xlog)-parameters.clipLimit*log(10)/20; % convert clipLimit from dB to log-magnitude
   Xlog(Xlog<Xfloor) = Xfloor;
end;

% Zero pad log spectrum:
LL = round(parameters.FFTsize/2); % samples
Xlog = [Xlog(1:LL);zeros((parameters.IFFTfactor-1)*parameters.FFTsize+1,1);Xlog(LL:-1:2)]; % ignore Nyquist frequency
Xlog = Xlog*parameters.IFFTfactor; % scale by IFFTfactor to ensure proper amplitude scaling in cepstral domain

% IFFT:
c = ifft(Xlog,'symmetric');

% Find cepstral peak in c in F0 range:
[CP,F0] = getPeak(c,fs,parameters.F0range,parameters.IFFTfactor);

return;

function w = createWindow(windowType,L)
% This function creates a window function vector of length L.

switch lower(windowType)
   case 'rectangular'
      w = ones(L,1);
   case 'hamming'
      w = hamming(L);
   case 'hanning'
      w = hanning(L);
   case 'blackman'
      w = blackman(L);
   case 'boxcar'
      w = ones(L,1);
end;

return;

function [CP,F0] = getPeak(c,fs,F0range,IFFTfactor)
% This function finds the peak in c at lags in the range of F0range.

% Account for IFFT factor in fs:
fs1 = fs*IFFTfactor; % sampling rate of c, after IFFT zero padding

% Convert F0range to index limits:
tRange = [ceil(fs1/F0range(2))+1,floor(fs1/F0range(1))+1]; % integer samples
tRange = [max(1,tRange(1)),min(length(c),tRange(2))]; % limit to indices into c

% Get peak of c in tRange:
[CP,CPindex] = max(c(tRange(1):tRange(2))); % index into [tRange(1):tRange(2)]
CPindex = CPindex + (tRange(1)-1); % index into c

% Get fundamental frequency at peak, account for IFFT expansion factor:
F0 = fs1/(CPindex-1); % convert index to lag with '-1'

return;

% Bye!
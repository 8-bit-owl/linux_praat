function [mLinear,mBark,X] = computeMoments(x,fs,parameters)
% This function calculates moments of audio data in x.  Moments are calculated from 4 spectral
% estimates: linear, bark, excitation, and loudness.
% Input:
%    x -- Nx1 real vector, WAV file data of speech
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters
%       .earphone -- ['er2'], used for calculation of specific loudness
%       .FFTsize -- integer scalar, FFT size [1024]
%       .windowType -- ['hamming'],'hanning','blackman' analysis window type
%       .windowSize -- real scalar, sec, analysis window size [20e-3]
%       .frameRate -- real scalar, frames/sec, analysis frames per second [100]
% Output:
%    mLinear, mBark -- struct of output variables
%       .m1-m4 -- 1xL real vectors, moments 1-4 for each frame, L frames
%       .S -- 1xL real vector, skewness of each frame
%       .K -- 1xL real vector, kurtosis of each frame
%       .f -- Kx1 real vector, frequencies of FFT bins of X, K bins
%    For mLinear, f in Hz, m1 in Hz, m2 in Hz^2, m3 in Hz^3, and m4 in Hz^4.  For mBark, Hz is
%    replaced with Bark in mLinear.
%    X -- KxL non-negative real matrix, K FFT bins in range [1,1+FFTsize/2], L frames in x,
%         magnitude squared power spectral density (not dB).

% Moments calculated according to Forrest et al. 1988.  Excitation curves and specific loudness
% curves from Moore et al. 1997. 
% References:
% Forrest, Weismer, Milenkovic, and Dougall, JASA vol. 84(1), pp. 115-123, July 1988
% Syrdal and Gopal, JASA vol. 79(4), pp. 1086-1100, April 1986 (Bark scale + low-freq corrections)
% Zwicker and Terhardt, JASA vol. 68(5), pp. 1523-1525, Nov. 1980 (Bark scale)
% Moore, Glasberg, and Baer, JAES, vol. 45(4), pp. 224-240, April 1997 (excitation/loudness)
%
% Mark Skowronski, October 16, 2012

% Parameters:
defaults(1).earphone = 'er2'; % earphone type used to calculate excitation/specific loudness curves
defaults(1).FFTsize = 1024;
defaults(1).windowType = 'hamming';
defaults(1).windowSize = 20e-3; % sec
defaults(1).frameRate = 100; % frames/sec

% Check inputs:
if nargin<3,
   parameters = defaults;
end;

% Check parameters for default entries:
parameters = checkParameters(parameters,defaults);

% Get power spectrum of frames of x:
X = getSpectrogram(x,fs,parameters); % magnitude squared power spectral density

% Set frequency of spectrogram bins:
f = [0:size(X,1)-1]'*(fs/parameters.FFTsize); % Hz, frequency of each FFT bin, linear frequency scale, COLUMN vector
fBark = getfBark(f); % converts Hz to Bark

% Calculate moments of X using f and fBark:
mLinear = calcMoments(X,f);
mBark = calcMoments(X,fBark);

% Get excitation and specific loudness curves:
%[excitationCurve,partialLoudnessCurve,freqERB] = processPartialLoudness(x,fs,earphone);

% Excitation moments:
%[E1,E2,E3,E4,ES,EK] = calcMoments(excitationCurve,freqERB); % excitationCurve in power units

% Specific loudness moments:
%[S1,S2,S3,S4,SS,SK] = calcMoments(partialLoudnessCurve,freqERB); % partialLoudnessCurve in sones/ERB (NOT power units)

return;


function output = calcMoments(P,f)
% This function calculates spectral moments from the spectrum (frequency weighted).
% Input:
%    P -- KxL non-negative real matrix, K FFT bins in range [1,1+FFTsize/2], L frames in x,
%         magnitude squared power spectral density.
%    f -- Kx1 real vector, frequency
%    P may be in power units/frequency or sones/frequency, frequency may be in Hz, Bark, or ERB.
% Output:
%    m1-m4 -- 1xL real vectors, moments 1-4 of P, scaled by f
%    S -- 1xL real vector, skewness of P
%    K -- 1xL real vector, kurtosis of P

% Weigth P by delta(f) (Forrest et al. 1988):
% Note: for bark scale, f saturates at 150 Hz, so delta(f)=0
deltaF = f(2:end)-f(1:end-1);
deltaF = [deltaF(1);deltaF]; % repeat deltaF(1) for first FFT bin
deltaF = repmat(deltaF,1,size(P,2)); % repeat for each column of P
P = P.*deltaF;

% Normalize each column of P:
p = P./repmat(sum(P,1),size(P,1),1); % p(k,l)>=0, sum(p,1)=1, p ~ probability density function

% Calculate first moment:
F = repmat(f,1,size(p,2)); % repeat for each frame
m1 = sum(F.*p,1);

% Calculate second-fourth moments:
M1 = repmat(m1,size(p,1),1); % repeat for each frequency bin
m2 = sum(((F-M1).^2).*p,1);
m3 = sum(((F-M1).^3).*p,1);
m4 = sum(((F-M1).^4).*p,1);

% Calculate skewness/kurtosis (unitless):
S = m3./(m2.^1.5);
K = m4./(m2.^2);

% Save output:
output(1).m1 = m1;
output(1).m2 = m2;
output(1).m3 = m3;
output(1).m4 = m4;
output(1).S = S;
output(1).K = K;
output(1).f = f;

return;


function fBark = getfBark(f)
% This function converts f in Hz to fBark in Bark.

% Apply low-frequency corrections to f (Syrdal & Gopal 1986, from Traunmuller 1981):
fhat = f;
fhat(f<150) = 150; % saturate at 150 Hz
fhat(f>=150 & f<200) = 0.8*f(f>=150 & f<200)+30; % slight compression in [150,200] Hz range
fhat(f>=200 & f<250) = 1.2*f(f>=200 & f<250)-50; % slight expansion in [200,250] Hz range

% Apply Bark scale transformation (Syrdal & Gopal 1986, from Zwicker & Terhardt 1980):
fBark = 13*atan(0.76*fhat/1000)+3.5*atan((fhat/7500).^2);

return;

function X = getSpectrogram(x,fs,parameters)
% This function generates a spectrogram.
% Input:
%    x -- Nx1 real vector, time domain signal
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of values that control spectrogram computation
% Output:
%    X -- KxL non-negative real matrix, K FFT bins in range [1,1+FFTsize/2], L frames in x,
%         magnitude squared power spectral density.

% Get window function:
w = getWindowFunction(fs,parameters);

% Divide x into frames:
xFrames = make_xFrames(x,fs,parameters);

% Apply window:
xWindowed = xFrames.*repmat(w,1,size(xFrames,2));

% Zero-mean each frame:
xWindowMean = mean(xWindowed,1); % mean of each COLUMN
xWindowed = xWindowed - repmat(xWindowMean,size(xWindowed,1),1);

% FFT windowed frames:
X = fft(xWindowed,parameters(1).FFTsize);
X = X(1:parameters(1).FFTsize/2+1,:); % lower half of frequencies
X = real(X).^2 + imag(X).^2; % magnitude squared

return;


function xFrames = make_xFrames(x,fs,parameters)
% This function breaks x up into a matrix of frames of data.

% Get length of window in samples:
L = getWindowLength(fs,parameters(1).windowSize); % samples per frame

% Determine number of frames, only frames with complete data are considered, no zero-padding last
% frame.  Note, frame counter used inside round() function to eliminate frame "drift" in cases of
% non-standard fs and/or frame rates.
numFrames = 1; % init with 1 frame
while L+round(numFrames*fs/parameters(1).frameRate) <= length(x), % check endpoint of (numFrames+1) frame
   numFrames = numFrames+1; % increment to next frame if endpoint <= length(x)
end;

% Init output:
xFrames = zeros(L,numFrames);

% Fill in output:
for p=1:numFrames,
   xTemp = x([1:L]+round((p-1)*fs/parameters(1).frameRate)); % may be ROW or COLUMN vector
   xFrames(1:L,p) = xTemp(:); % COLUMN vector
end;

return;


function L = getWindowLength(fs,windowSize)
% This function converts windowSize from sec to samples.

L = round(fs*windowSize); % samples

return;


function w = getWindowFunction(fs,parameters)
% This function returns a COLUMN vector of a window function useful for spectral analysis.

L = getWindowLength(fs,parameters(1).windowSize);
switch lower(parameters(1).windowType)
   case 'blackman'
      w = blackman(L);
   case 'hanning'
      w = hanning(L);
   otherwise % Hamming window is default
      w = hamming(L);
end;

return;



% Bye!

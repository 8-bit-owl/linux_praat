function [output,parameters,X,FBout] = HFCC(x,fs,parameters)
% This function implements the algorithm for extracting human factor cepstral coefficients (HFCCs),
% similar to Skowronski and Harris, JASA, 2004.  Filter bank design is simplified.
% Input [default values]:
%  x -- Nx1 real vector, time domain signal from which HFCCs are extracted
%  fs -- real scalar, Hz, sampling rate of x [16000]
%  parameters -- struct of feature extraction parameters,
%     .fb -- numFilters x FFTsize/2+1 real matrix, HFCC filter bank matrix (constructed if empty/missing)
%     .DCTmatrix -- numCoeffs x numFilters real matrix, DCT matrix (constructed if empty/missing)
%     .FFTsize -- integer scalar, length of FFT used to estimate spectrum of each frame [1024]
%     .numFilters -- integer scalar, number of HFCC filters to use in filter bank [30]
%     .numCoeffs -- integer scalar, number of coefficients returned for each frame [13]
%     .freqRange -- 1x2 real vector, Hz, [min max] center frequencies of HFCC filter bank [70 7000]
%     .Efactor -- real scalar, ERB scale factor, affects bandwidth of HFCC filters [1.0]
%     .deltaSize -- integer scalar, +/- number of frames to use to estimate local slope [1]
%     .windowType -- string of FFT window name: ['hamming'],'hanning','blackman'
%     .applyPreemphasis -- logical scalar: 1=apply pre-emphasis; 0=don't [1]
%     .preemphasisAlpha -- real scalar, FIR zero location of pre-emphasis filter [0.95]
%     .windowSize -- real scalar, sec, duration of analysis window [20e-3]
%     .frameRate -- real scalar, frames/sec, number of analysis windows per second [100]
%     .applyCMS -- logical scalar, 1=apply cepstral mean subtraction; 0=don't [1]
%     .fs -- real scalar, Hz, sampling rate used to generate fb (parameters.fb redesigned if fs ~=
%            parameters.fs) [16000]
%     .ccMethod -- string, ['HFCC'],'MFCC' method for calculating cepstral coefficients
% Output:
%  output -- struct of outputs
%     .cc -- numCoeffs x numFrames real matrix, cepstral coefficients
%     .Dcc -- numCoeffs x numFrames real matrix, delta cepstral coefficients
%     .DDcc -- numCoeffs x numFrames real matrix, delta delta cepstral coefficients
%     .time -- 1xnumFrames real vector, sec, time at start of each frame
%  parameters -- same as input variable, with fb and DCTmatrix constructed if missing from input or
%                fs mismatch
%  X -- N x numFrames real matrix, spectrogram magnitude of x
%  FBout -- numFilters x numFrames real matrix, filter bank output, log scaled

% Mark Skowronski, October 4, 2011
% Constructed from fbInit.m (2003) and getFeatures.m (2002).

% Create defaults:
defaults = struct([]); % add fb and DCTmatrix below
defaults(1).FFTsize = 1024;
defaults(1).numFilters = 30;
defaults(1).numCoeffs = 13;
defaults(1).freqRange = [70 7000]; % Hz
defaults(1).Efactor = 1.0;
defaults(1).deltaSize = 1;
defaults(1).windowType = 'hamming';
defaults(1).applyPreemphasis = true;
defaults(1).preemphasisAlpha = 0.95;
defaults(1).windowSize = 20e-3; % sec
defaults(1).frameRate = 100; % frames/sec
defaults(1).applyCMS = true;
defaults(1).fs = 16000; % Hz
defaults(1).ccMethod = 'HFCC';

% Check inputs:
if nargin<2,
   fs = 16000; % Hz
end;
if nargin<3,
   parameters = defaults;
end;

% Fill in missing parameter fields with defaults:
parameters = checkParameters(parameters,defaults);

% Make filters if necessary:
parameters = checkFilters(fs,parameters);

% Apply pre-emphasis if necessary:
if logical(parameters(1).applyPreemphasis),
   x = applyPreemphasis(x,parameters(1).preemphasisAlpha);
end;

% Construct spectrogram for x:
[X,tFrames] = getSpectrogram(x,fs,parameters);

% Get log power from X:
logPower = log10(mean(X.^2,1)); % mean along each COLUMN.

% Combine X with fb and DCTmatrix to get cepstral coeffs:
DCTmatrixShort = parameters(1).DCTmatrix(2:end,:); % ignore first row (DC) of DCT matrix
FBout = log10(parameters(1).fb*X); % filter bank output, log scaled
cc = DCTmatrixShort*FBout; % numCoeffs x numFrames matrix

% Prepend log power to cc:
cc = [logPower;cc];

% Apply CMS if necessary:
if logical(parameters(1).applyCMS),
   cc = applyCMS(cc);
end;

% Get Dcc, DDcc:
Dcc = getDeltaCoeffs(cc,parameters(1).deltaSize);
DDcc = getDeltaCoeffs(Dcc,parameters(1).deltaSize);

% Save:
output(1).cc = cc;
output(1).Dcc = Dcc;
output(1).DDcc = DDcc;
output(1).time = tFrames;

return;

function cc = applyCMS(cc)
% This function applies cepstral mean subtraction to the matrix of cepstral coeffs.

% Calculate mean of each ROW:
ccMean = mean(cc,2); % COLUMN vector

% Subtract mean from cc:
cc = cc - repmat(ccMean,1,size(cc,2));

return;

function [X,tFrames] = getSpectrogram(x,fs,parameters)
% This function generates a spectrogram useful for HFCC analysis.

% Get window function:
w = getWindowFunction(fs,parameters);

% Divide x into frames:
[xFrames,tFrames] = make_xFrames(x,fs,parameters);

% Apply window:
xWindowed = xFrames.*repmat(w,1,size(xFrames,2));

% Zero-mean each frame:
xWindowMean = mean(xWindowed,1); % mean of each COLUMN
xWindowed = xWindowed - repmat(xWindowMean,size(xWindowed,1),1);

% FFT windowed frames:
X = abs(fft(xWindowed,parameters(1).FFTsize)); % magnitude

% Return only lower half of frequencies:
X = X(1:parameters(1).FFTsize/2+1,:);

return;

function [xFrames,tFrames] = make_xFrames(x,fs,parameters)
% This function breaks x up into a matrix of frames of data and also returns the time of the
% beginning of each frame.

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
tFrames = zeros(1,numFrames);
% Fill in output:
for p=1:numFrames,
   tTemp = [1:L]+round((p-1)*fs/parameters(1).frameRate); % index into x
   tFrames(p) = tTemp(1); % index into x
   xTemp = x(tTemp); % may be ROW or COLUMN vector
   xFrames(1:L,p) = xTemp(:); % COLUMN vector
end;

% Convert tFrames from index to sec:
tFrames = (tFrames-1)/fs;

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

function parameters = checkFilters(fs,parameters)
% This function checks to see if the HFCC filter bank and DCT matrix exist and creates them if
% needed.

makeFB = true; % make fb and DCTmatrix unless present
if all(isfield(parameters,{'fb','DCTmatrix'})), % if fields are present
   if ~isempty(parameters(1).fb) && ~isempty(parameters(1).DCTmatrix), % fields are not empty
      if parameters.fs == fs, % input fs and design fs are identical
         makeFB = false; % so don't fb or DCTmatrix
      end;
   end;
end;
if makeFB,
   parameters = makeHFCCfilters(parameters,fs);
end;

return;

function x = applyPreemphasis(x,alphaZero)
% This function applies a pre-emphasis HPF to x, first-order FIR with zero at z = alphaZero.

% Make filter:
h = [1 -alphaZero];

% Filter:
x = filter(h,1,x);

return;

function parameters = makeHFCCfilters(parameters,fs)
% This function constructs the HFCC filter bank and DCT matrix from parameters and fs.  The original
% filter bank design set the center frequencies of the first/last filters such that the filters
% taper to the min/max frequencies in the specified freqRange parameter.  This constraint
% complicates the calculation of center freqs.  Furthermore, center freqs change as Efactor changes,
% and center freqs converge to a single value when Efactor is large enough.  Also, filter endpoints
% are equally spaced around the center frequency in mel frequency; when warped to linear frequency,
% frequency endpoints are such that ERB = (fHigh-fLow)/2.  Triangle legs are straight in linear
% frequency.
%
% To simplify center frequency calculations, the min/max frequencies in freqRange specify the
% frequencies of the first/last filter.  Thus, the first/last filters taper to frequencies beyond
% the specified min/max freqs in freqRange.  Filters are truncated at 0 Hz and fs/2 Hz.  Also,
% filter endpoints are equally spaced around the center frequency in linear frequency; thus, fLow =
% fCenter - ERB and fHigh = fCenter + ERB such that ERB = (fHigh-fLow)/2 as before.
%
% The function also supports the option to create the Davis and Mermelstein filterbank: first 10
% filter center frequencies = [100:100:1000] Hz, higher center frequencies are spaced 5 per octave.
% Flower of first filter = 0 Hz, Fupper of the last filter is the largest Fupper below the Nyquist
% rate (or 8 kHz if Nyquist rate > 8 kHz).

% Design center frequency and lower/upper frequency for each filter:
switch lower(parameters.ccMethod),
   case 'hfcc'
      % Convert freqRange from linear frequency to mel frequency:
      freqRangeMel = lin2mel(parameters.freqRange);

      % Calculate filter center frequencies in mel frequency:
      N = parameters.numFilters;
      fCenterMel = freqRangeMel(1) + [0:N-1]/(N-1)*(freqRangeMel(2)-freqRangeMel(1));

      % Convert center frequencies to linear frequency:
      fCenter = mel2lin(fCenterMel); % Hz

      % Get ERB bandwidths for center frequencies:
      BW = getERB(fCenter); % Hz

      % Scale bandwidth by ERB scale factor:
      BW = BW*parameters.Efactor;
      
      % Set fLower/fUpper:
      fLower = fCenter - BW;
      fUpper = fCenter + BW;
      numFilters = parameters.numFilters;
   case 'mfcc'
      % Determine number of center frequencies (24 filters for fs >= 16000 Hz):
      % Note: 1000*(2^1/5)^p <= min(8000,fs/2), p upper frequencies above 1 kHz that are <= 8 kHz or Nyquist rate.
      numFilters = 10 + floor(log(min(8000,fs/2)/1000)/(1/5*log(2)))-1; % 10 between 100-1000 Hz, last center frequency = p-1
      
      % Design center frequencies:
      fCenter = zeros(1,numFilters);
      fCenter(1:10) = [100:100:1000]; % linear spacing
      fCenter(11:end) = 1000*((2^(1/5)).^[1:numFilters-10]);
      
      % Create fLower/fUpper:
      fLower = [0,fCenter(1:end-1)];
      fUpper = [fCenter(2:end),fCenter(end)*(2^(1/5))];
end;

% Design filter bank matrix:
parameters.fb = getHFCCfb(fCenter,fLower,fUpper,fs,parameters.FFTsize);

% Design DCT matrix:
parameters.DCTmatrix = getDCTmatrix(parameters.numCoeffs,numFilters);

% Save fs used to create fb:
parameters(1).fs = fs;

return;

function DCTmatrix = getDCTmatrix(numCoeffs,numFilters)
% This function constructs a type II DCT matrix:
% X(k) = sqrt(2/N)*sum_n=0^N-1 x(n)*cos(pi/N*(n+1/2)*k), N=numFilters, k=0:numCoeffs-1
% In HFCC, X(0) is replaced with log energy, but the DCT matrix returned still contains X(0) in the
% first row for completeness.
% Input:
%  numCoeffs -- integer scalar, number of cepstral coefficients output from matrix operation
%  numFilters -- integer scalar, number of filters in HFCC filter bank
% Output:
%  DCTmatrix -- numCoeffs x numFilters real matrix, orthogonal: DCTmatrix*DCTmatrix' = I

% Make time and frequency vectors:
n = [0:numFilters-1]; % ROW vector, time
k = [0:numCoeffs-1]'; % COLUMN vector, frequency

% Make DCT matrix:
DCTmatrix = cos(pi/numFilters*k*(n+1/2)); % numCoeffs x numFilters matrix

% Make orthogonal:
DCTmatrix = DCTmatrix * sqrt(2/numFilters);
DCTmatrix(1,:) = DCTmatrix(1,:)/sqrt(2); % DC row needs extra scaling

return;

function fb = getHFCCfb(fCenter,fLower,fUpper,fs,FFTsize)
% This function creates the HFCC filter bank.
% Input:
%  fCenter -- 1xN real vector, Hz, center frequency of N filters in bank, ordered min to max
%  BW -- 1xN real vector, Hz, bandwidth of filters at fCenter
%  fs -- real scalar, Hz, sampling rate
%  FFTsize -- integer scalar, total length of FFT
% Output:
%  fb -- NxL real matrix, filter bank coefficients, L = FFTsize/2+1, where bin L of FFT corresponds
%  to fs/2.

% Init output:
N = length(fCenter);
L = FFTsize/2+1;
fb = zeros(N,L);

% Create each filter:
for p=1:N,
   fb(p,:) = makeHFCCfilter(fCenter(p),fLower(p),fUpper(p),fs,L);
end;

return;

function fbRow = makeHFCCfilter(fCenter,fLow,fHigh,fs,L)
% This function makes a single HFCC filter (one row in fb).

% Init output:
fbRow = zeros(1,L); % filter gain for FFT bins 1 to L

% Find ranges of FFT bins in lower/upper triangle halves:
f = [0:L-1]/(L-1)*fs/2; % Hz, frequencies of corresponding FFT bins 1 (0 Hz) to L (fs/2 Hz)
fLowRange = f>=fLow & f<=fCenter; % 1 for FFT bins in lower half of triangle filter, 0 otherwise
fHighRange = f>=fCenter & f<=fHigh; % 1 for FFT bins in upper half of triangle filter, 0 otherwise

% Create lower/upper triangle legs:
fLowLeg = (f-fLow)/(fCenter-fLow); % 0 at f=fLow, 1 at f=fCenter
fHighLeg = (f-fHigh)/(fCenter-fHigh); % 0 at f=fHigh, 1 at f=fCenter

% Save lower/upper legs in fbRow (zero, otherwise):
fbRow(fLowRange) = fLowLeg(fLowRange);
fbRow(fHighRange) = fHighLeg(fHighRange);

return;

function BW = getERB(f)
% This function calculates equivalent rectangular bandwidth for linear frequencies in f (Hz)
% according to Moore and Glasberg, JASA, 1983.

BW = 6.23*(f/1000).^2 + 93.39*(f/1000) + 28.52; % Hz

return;

function fMel = lin2mel(f)
% This function uses O'Shaughnessy's equation to convert linear frequency to mel frequency.

A = 1000/log10(1+1000/700); % ~2595, scale factor such that 1000 Hz --> 1000 mel
fMel = A*log10(1+f/700);

return;

function f = mel2lin(fMel)
% This function is the inverse of lin2mel.

A = 1000/log10(1+1000/700); % ~2595, scale factor such that 1000 Hz --> 1000 mel
f = (10.^(fMel/A)-1)*700;

return;

% Bye!
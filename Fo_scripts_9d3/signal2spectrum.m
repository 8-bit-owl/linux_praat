function [X,f,t,parameters] = signal2spectrum(x,fs,parameters)
% This function converts x into a frame-based spectrum of one of the following types:
% -- power spectral density, linear frequency scale
% -- power spectral density, bark frequency scale
% -- power spectral density, ERBrate frequency scale
% -- excitation pattern, equally spaced center frequencies in ERBrate space
% -- specific loudness pattern, equally spaced center frequencies in ERBrate space
% Input:
%    x -- Jx1 real vector, time domain signal to analyze, J samples
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of relevant parameters
%       .spectrumType -- ['PSDlinear'],'PSDbark','PSDERBrate','Excitation','SpecificLoudness'
%       .transducer -- string, name of acoustic transducer, ['freefield'],'er2', see
%          outerEarConversion.m, applied to excitation and specific loudness only
%       .FFTsize -- integer scalar, length of FFT used to estimate spectrum of each frame [1024]
%       .windowType -- string of FFT window name: ['hamming'],'hanning','blackman','rectangular'
%       .windowSize -- real scalar, sec, duration of analysis window [20e-3]
%       .frameRate -- real scalar, frames/sec, number of analysis windows per second [100]
%       .fcRange -- 1x2 real vector, Hz, [min max] frequencies of excitation filter bank center
%          frequencies [50 7500]
%       .ERBspacing -- real scalar, ERBrate units, spacing between adjacent center frequencies
%          (equal spacing in ERBrate frequency space) [0.1]
%       .auditoryFilterType -- string, ['GCC'],'roex', name of auditory filter method to use to
%          calculate excitation pattern
%       .zeroMean -- [false],true logical flag, apply zero mean to each frame before spectrum
%          calculation
% Output:
%    X -- MxT real matrix, spectrum matrix, M frequency bins, T analysis frames, each spectrum in
%         one COLUMN of X
%    f -- Mx1 real vector, frequency units depending on spectrumType
%    t -- 1xT real vector, sec, T analysis frames, time of BEGINNING of each analysis frame
%    parameters -- same as input, or defaults if empty. Includes FFTsizeOutput which may be larger
%       than FFTsize in case FFTsize < analysis window length.
%
% Note: to analyze x in a single analysis frame (ignoring windowSize), set frameRate = 0 or set
% windowSize = length(x) and frameRate <= fs.
% Note: make sure to scale x to desired RMS value (dB SPL) for excitation and specific loudness
% calculation, as excitation filter bandwidth is level-dependent and specific loudness calculations
% are level-dependent.

% Mark Skowronski, May 20, 2013

% Set defaults:
defaults = struct([]); % add fb and DCTmatrix below
defaults(1).spectrumType = 'PSDlinear';
defaults(1).transducer = 'freefield';
defaults(1).FFTsize = 1024;
defaults(1).windowType = 'hamming';
defaults(1).windowSize = 20e-3; % sec
defaults(1).frameRate = 100; % frames/sec
defaults(1).fcRange = [50 7500]; % Hz
defaults(1).ERBspacing = 0.1; % ERBrate units
defaults(1).auditoryFilterType = 'GCC'; % compressed gammachirp auditory filterbank
defaults(1).zeroMean = false;

% Init output:
X = [];
f = [];
t = [];

% Check input:
if nargin<3,
   parameters = defaults;
else % set default fields if necessary
   parameters = checkParameters(parameters,defaults);
end;
if nargin<1,
   return; % just return parameters and empty output variables
end;

% Separate x into frames, convert t from samples to sec:
[xFrames,tSample] = makeFrames(x,fs,parameters); % cell array
t = (tSample-1)/fs; % sec

% Zero mean each frame if necessary:
if parameters.zeroMean
   for p=1:length(xFrames),
      xFrames{p} = xFrames{p}-mean(xFrames{p});
   end;
end;

% Construct spectrum from frames of x:
switch lower(parameters.spectrumType)
   case 'psdlinear',
      [X,f] = makePSDMatrix(xFrames,fs,parameters); % f in Hz, FFT bin frequencies
   case 'psdbark',
      [X,f] = makePSDMatrix(xFrames,fs,parameters);
      f = f2bark(f); % Bark
   case 'psderbrate',
      [X,f] = makePSDMatrix(xFrames,fs,parameters);
      f = f2ERBrate(f); % ERBrate
   case 'excitation',
      [X,f] = makeExcitationMatrix(xFrames,fs,parameters); % f in Hz, filterbank center frequencies
   case 'specificloudness',
      [X,f] = makeSpecificLoudnessMatrix(xFrames,fs,parameters); % f in Hz, filterbank center frequencies
end;

return;

function [X,f] = makeSpecificLoudnessMatrix(xFrames,fs,parameters)
% This function calculates the specific loudness pattern for each frame of x.
% Input:
%    xFrames -- Tx1 cell array, frame of data in each cell, T frames
%    fs -- real scalar, Hz, sampling rate of data in each cell
%    parameters -- struct of parameters
%       .fcRange -- 1x2 real vector, Hz, [min max] frequencies of center frequency range
%       .ERBspacing -- real scalar, ERB/ERB, spacing between adjacent center frequencies in ERBrate
%       .FFTsize -- FFT length
%       .windowType -- string of window type
%       .auditoryFilterType -- string, ['GCC'],'roex' of auditory filter method to use
%       .transducer -- string, name of acoustic transducer, ['freefield'],'er2', see
%                      outerEarConversion.m
% Output:
%    X -- QxT real matrix, power units, excitation of analysis frames, Q filterbank center
%         frequencies
%    f -- Qx1 real vector, Hz, filterbank center frequencies
% Note: if FFTsize is less than the number of samples in a frame L, then FFTsize = L in getSpectrum
% function.
% Note: a filter bank is designed for each analysis frame (slow).

% Create excitation pattern for each frame, center frequencies:
[X1,f] = makeExcitationMatrix(xFrames,fs,parameters); % f in Hz, filterbank center frequencies

% Create specific loudness pattern for each excitation pattern:
X = zeros(size(X1)); % init output
Enoise = 0; % assume zero masking noise of signal
for p=1:size(X1,2), % for each analysis frame
   [junk,junk2,X(:,p)] = excitation2loudness(X1(:,p),Enoise,f); % ignore signal and noise total loudness values
end;

return;


function [X,f] = makeExcitationMatrix(xFrames,fs,parameters)
% This function calculates the excitation pattern for each frame of x.
% Input:
%    xFrames -- Tx1 cell array, frame of data in each cell, T frames
%    fs -- real scalar, Hz, sampling rate of data in each cell
%    parameters -- struct of parameters
%       .fcRange -- 1x2 real vector, Hz, [min max] frequencies of center frequency range
%       .ERBspacing -- real scalar, ERB/ERB, spacing between adjacent center frequencies in ERBrate
%       .FFTsize -- FFT length
%       .windowType -- string of window type
%       .auditoryFilterType -- string, ['GCC'],'roex' of auditory filter method to use
%       .transducer -- string, name of acoustic transducer, ['freefield'],'er2', see
%                      outerEarConversion.m
% Output:
%    X -- QxT real matrix, power units, excitation of analysis frames, Q filterbank center
%         frequencies
%    f -- Qx1 real vector, Hz, filterbank center frequencies
% Note: if FFTsize is less than the number of samples in a frame L, then FFTsize = L in getSpectrum
% function.
% Note: a filter bank is designed for each analysis frame (slow).

% Make filterbank center frequencies:
f = makeFc(parameters.fcRange,parameters.ERBspacing);

% Get first excitation, init output:
[X1,junk,parameters] = signal2excitation(xFrames{1},fs,f,parameters,[],true); % check parameters
X = zeros(length(X1),length(xFrames));
X(:,1) = X1;

% Process all other frames:
for p=2:length(xFrames),
   X(:,p) = signal2excitation(xFrames{p},fs,f,parameters,[],false); % design FB for each frame
end;

return;

function [X,f] = makePSDMatrix(xFrames,fs,parameters)
% This function calculates the power spectral density for each frame of x.
% Input:
%    xFrames -- Tx1 cell array, frame of data in each cell, T frames
%    fs -- real scalar, Hz, sampling rate of data in each cell
%    parameters -- struct of parameters
%       .FFTsize -- FFT length
%       .windowType -- string of window type
% Output:
%    X -- QxT real matrix, power units, PSD of analysis frames, Q FFT bins (lower half of spectrum)
%    f -- Qx1 real vector, Hz, frequencies of FFT bins
% Note: if FFTsize is less than the number of samples in a frame L, then FFTsize = L in getSpectrum
% function.

% Get first spectrum and f, init output:
[X1,f] = getSpectrum(xFrames{1},fs,parameters.FFTsize,parameters.windowType);
X = zeros(length(X1),length(xFrames));
X(:,1) = X1;

% Process all other frames:
for p=2:length(xFrames),
   s = getSpectrum(xFrames{p},fs,parameters.FFTsize,parameters.windowType);
   X(1:length(s),p) = s; % may be shorter than X1 if FFTsize too small
end;

return;


% Bye!
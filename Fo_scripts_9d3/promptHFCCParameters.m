function parameters = promptHFCCParameters(parameters,defaultsOnly)
% This function prompts for HFCC parameters or returns default values according to input flag.
% Input:
%    parameters -- scalar struct or empty, any input fields
%    defaultsOnly -- logical scalar, [false]: prompt user for parameter values, true: skip prompt
%       and return default values
% Output:
%    parameters -- scalar struct, same as input with added fields
%       .ccMethod -- string, ['HFCC'],'MFCC'
%       .windowType -- string, ['Hamming'],'Hanning','Blackman','Rectangular' for CC spectral
%          estimation 
%       .windowSize -- real scalar, sec, analysis window length [0.02]
%       .frameRate -- real scalar, fps, analysis frames per second [100]
%       .FFTsize -- integer scalar, samples, for spectral estimation [1024]
%       .numFilters -- integer scalar, number of triangular filters in filter bank [40]
%       .numCoeffs -- integer scalar, number of cepstral coefficients estimated from each frame [13]
%       .freqRange -- real 1x2 vector, Hz, [min max] frequencies of filter bank [70 7500]
%       .Efactor -- real scalar, triangular filter bandwidth scale factor [1]
%       .deltaSize -- integer scalar, frames, +/- number of frames used to estimate delta
%          coefficients [1]
%       .applyPreemphasis -- logical scalar, [true],false for applying pre-emphasis
%       .preemphasisAlpha -- real scalar, pre-emphasis filter zero location [0.95]
%       .applyCMS -- logical scalar, [true],false for applying cepstral mean subtraction
%       .fsHFCC -- real scalar, Hz, RESAMPLED sampling rate of x, ensures fair comparison of
%          measures from data files with various original sampling rates [16000]
%       .HPF_order -- integer scalar, order of high-pass filter applied to x after resampling, set
%          to 0 to skip high-pass filter operation [10] 
%       .HPF_fc -- real scalar, Hz, cutoff frequency of high-pass filter applied after resampling,
%          set to 0 to skip high-pass filter operation [70] 
%       .trimThreshold -- positive real scalar, dB, beginning/end of x energy contour below
%          threshold trimmed off [25]
%       .trimTimeConstant -- positive real scalar, sec, energy contour time constant (exponential
%          moving average estimate), set to 0 to skip endpoint trimming [0.02]
%       .calcASL -- logical scalar, true,[false] for calculating active speech level values
%       .calcSWIPE -- logical scalar, true,[false] for calculating SWIPE pitch and strength

% Set defaults:
defaults(1).ccMethod = 'HFCC';
defaults(1).windowType = 'Hamming';
defaults(1).windowSize = 20e-3; % sec
defaults(1).frameRate = 100; % frames/sec
defaults(1).FFTsize = 1024;
defaults(1).numFilters = 40;
defaults(1).numCoeffs = 13;
defaults(1).freqRange = [70 7500]; % Hz
defaults(1).Efactor = 1;
defaults(1).deltaSize = 1;
defaults(1).applyPreemphasis = true;
defaults(1).preemphasisAlpha = 0.95;
defaults(1).applyCMS = true;
defaults(1).fsHFCC = 16000;
defaults(1).HPF_order = 10; % HPF applied in processSentence.m
defaults(1).HPF_fc = 70; % Hz
defaults(1).trimThreshold = 25; % dB, endpoint-trimming applied in processSentence.m > trimEndpoints.m
defaults(1).trimTimeConstant = 20e-3; % sec
defaults(1).calcASL = false; % don't measure active speech level (slow)
defaults(1).calcSWIPE = false; % don't measure SWIPE pitch (slow)

% Check inputs:
if nargin<2,
   defaultsOnly = false; % prompt user for inputs
end;
if nargin<1 || isempty(parameters), % set defaults
   parameters = defaults;
else % Ensure all default fields exist in parameters:
   parameters = checkParameters(parameters,defaults);
end;

% Set display strings for special cases:
applyPreemphasisDisplay = 'false';
if parameters.applyPreemphasis,
   applyPreemphasisDisplay = 'true';
end;
applyCMSDisplay = 'false';
if parameters.applyCMS,
   applyCMSDisplay = 'true';
end;
calcASLDisplay = 'false';
if parameters.calcASL,
   calcASLDisplay = 'true';
end;
calcSWIPEDisplay = 'false';
if parameters.calcSWIPE,
   calcSWIPEDisplay = 'true';
end;

% Set parameter prompts and default values:
parameterBlock = {...
   'CC method [HFCC, MFCC]',parameters.ccMethod;...
   'Window type [Hamming, Hanning, Blackman, Rectangular]',parameters.windowType;...
   'Window size, ms',num2str(parameters.windowSize*1e3);...
   'Frame rate, fps',num2str(parameters.frameRate);...
   'FFT size, samples',num2str(parameters.FFTsize);...
   'Number filters',num2str(parameters.numFilters);...
   'Number coefficients',num2str(parameters.numCoeffs);...
   'MINIMUM frequency range, Hz',num2str(parameters.freqRange(1));...
   'MAXIMUM frequency range, Hz',num2str(parameters.freqRange(2));...
   'E factor',num2str(parameters.Efactor);...
   'Delta size, frames',num2str(parameters.deltaSize);...
   'Apply pre-emphasis [true, false]',applyPreemphasisDisplay;...
   'Pre-emphasis alpha',num2str(parameters.preemphasisAlpha);...
   'Apply cepstral mean subtraction [true, false]',applyCMSDisplay;...
   'RESAMPLED sampling rate, Hz [data resampled from original sampling rate]',num2str(parameters.fsHFCC);...
   'High-pass filter order [0 to disable filter]',num2str(parameters.HPF_order);...
   'High-pass filter cut-off frequency, Hz',num2str(parameters.HPF_fc);...
   'Endpoint trim threshold, dB',num2str(parameters.trimThreshold);...
   'Endpoint trim time constant, ms [0 to skip endpoint trimming]',num2str(parameters.trimTimeConstant*1e3);...
   'Calculate active speech level [true, false]',calcASLDisplay;...
   'Calculate SWIPE during HFCC [true, false]',calcSWIPEDisplay};

% Split into smaller groups:
displayString1 = parameterBlock(1:10,1);
displayString2 = parameterBlock(11:end,1);
valueString1 = parameterBlock(1:10,2);
valueString2 = parameterBlock(11:end,2);

% Prompt for parameter settings if necessary:
if ~defaultsOnly, % prompt user for values
   a1 = inputdlg(displayString1,'CC PARAMETERS 1/2',[1 80],valueString1); % [] if cancel button
   a2 = inputdlg(displayString2,'CC PARAMETERS 2/2',[1 80],valueString2); % [] if cancel button
else % use default values without prompting
   a1 = []; % same as cancel button
   a2 = [];
end;

% Handle cancel button:
if isempty(a1), % cancel button
   a1 = valueString1; % use initial values
end;
if isempty(a2), % cancel button
   a2 = valueString2; % use initial values
end;
a = [a1;a2];

% Store input parameters:
parameters(1).ccMethod = a{1};
parameters(1).windowType = a{2};
parameters(1).windowSize = str2num(a{3})/1000; % sec
parameters(1).frameRate = str2num(a{4}); % frames/sec
parameters(1).FFTsize = str2num(a{5});
parameters(1).numFilters = str2num(a{6});
parameters(1).numCoeffs = str2num(a{7});
parameters(1).freqRange = [str2num(a{8}) str2num(a{9})]; % Hz
parameters(1).Efactor = str2num(a{10});
parameters(1).deltaSize = str2num(a{11});
parameters(1).applyPreemphasis = str2num(a{12});
parameters(1).preemphasisAlpha = str2num(a{13});
parameters(1).applyCMS = str2num(a{14});
parameters(1).fsHFCC = str2num(a{15});
parameters(1).HPF_order = str2num(a{16}); % HPF applied in processSentence.m
parameters(1).HPF_fc = str2num(a{17}); % Hz
parameters(1).trimThreshold = str2num(a{18}); % dB, endpoint-trimming applied in processSentence.m > trimEndpoints.m
parameters(1).trimTimeConstant = str2num(a{19})/1000; % sec
parameters(1).calcASL = str2num(a{20}); % don't measure active speech level (slow)
parameters(1).calcSWIPE = str2num(a{21}); % don't measure SWIPE during HFCC (slow)

return;
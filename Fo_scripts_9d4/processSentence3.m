function [output,parameters,x,X] = processSentence3(x,fs,parameters)
% This function processes an utterance of a sentence for HFCC variation measures.  The following
% steps are followed:
% 1) Resample x to 16 kHz for the sake of cepstral analysis,
% 2) High-pass filter at 70 Hz to remove low-frequency noise,
% 3) Trim endpoints (silence regions) at beginning/end of signal using energy envelope,
% 4) Calculate active speech level and activity factor (ITU-T P.56 Method B),
% 5) Calculate HFCCs and delta coefficients,
% 6) Calculate HFCC standard deviation sum.
%
% Input:
%    x -- Nx1 real vector, time series of utterance
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters used to process x [default values]
%    Note: see HFCC.m for HFCC parameters.
%       .fsHFCC -- real scalar, Hz, resampled rate of x [16000]
%       .HPF_order -- integer scalar, order of high-pass filter for removing low-freq. noise [10]
%       .HPF_fc -- real scalar, Hz, cutoff frequency of high-pass filter
%       .trimThreshold -- real scalar, dB, energy threshold for endpoint trimming [25]
%       .trimTimeConstant -- real scalar, sec, exp. moving average time constant [20e-3]
%       .calcASL -- logical scalar, true = calculate active speech level, false = don't [true]
%    Note: see HFCC.m for HFCC parameters.  Set parameters.fb to [] to design filter bank and DCT
%          matrix.
% Output:
%    output -- struct of output variables
%       .fs -- real scalar, Hz, resampled sampling rate, same as fsHFCC
%       .lenXtrimmed -- real scalar, sec, trimmed length of x
%       .ASL -- struct of active speech level function output (see ASL.m)
%       .activeSpeechLevel -- real scalar, dB DFS, active speech level of x
%       .activityFactor -- real scalar, unitless, ratio of active duration and total duration of x
%       .HFCC -- struct of HFCC output (see HFCC.m)
%       .ccStdSum -- real scalar, unitless, HFCC standard deviation sum measure
%       .DccStdSum -- real scalar, unitless, delta HFCC standard deviation sum measure
%       .DDccStdSum -- real scalar, unitless, delta delta HFCC standard deviation sum measure
%       .euclDistMean -- real scalar, unitless, average Euclidean distance between adjacent
%                          frames
%       .absDistMean -- real scalar, unitless, average absolute distance between adjacent frames
%       .DccMeasures -- same as ccMeasures, but for delta HFCC features
%       .DDccMeasures -- same as ccMeasures, but for delta-delta HFCC features
%       .FBMeasures, .DFBMeasures, .DDFBMeasures -- same as CC measures, calculated from filter bank
%          output
%    parameters -- same as input, but with added default fields
%    x -- Mx1 real vector, trimmed version of input x, M<=N
%    X -- spectrogram of x

% Mark Skowronski, October 2, 2012
% Edited by Lisa Kopf on September 14, 2012

% Set defaults for sentence preparations:
defaults(1).fsHFCC = 16000; % Hz
defaults(1).HPF_order = 10;
defaults(1).HPF_fc = 70; % Hz
defaults(1).trimThreshold = 25; % dB
defaults(1).trimTimeConstant = 20e-3; % sec
defaults(1).calcASL = true; 

% Check inputs:
if nargin<3,
   parameters = defaults;
else % add defaults fields to parameters if absent
   parameters = checkParameters(parameters,defaults);
end;

% Resample to fsHFCC:
[num,den] = rat(parameters.fsHFCC/fs);
x = resample(x,num,den);
output(1).fs = parameters.fsHFCC; % Hz
fs = parameters.fsHFCC;

% HPF x using 2nd-order filtering:
[zHPF,pHPF,kHPF] = butter(parameters.HPF_order,parameters.HPF_fc/(fs/2),'high'); % returns poles/zeros/gain
x = xFilter(x,zHPF,pHPF,kHPF);

% Trim endpoints of x:
x = trimEndpoints(x,fs,parameters.trimThreshold,parameters.trimTimeConstant);
output(1).lenXtrimmed = length(x)/fs; % sec

% Get active speech level:
output(1).ASL = ASL(x,fs); % use default parameters
output(1).activeSpeechLevel = output(1).ASL.activeSpeechLevel;
output(1).activityFactor = output(1).ASL.activityFactor;

% % Get HFCCs:
% [output(1).HFCC,parameters,X] = HFCC(x,fs,parameters);

% Get HFCCs and filter bank output:
[output(1).HFCC,parameters,X,FBout] = HFCC(x,fs,parameters); % Designs filterbank and DCT matrix on first call

% Calculate delta and delta-delta for FBout:
DFBout = getDeltaCoeffs(FBout,parameters.deltaSize);
DDFBout = getDeltaCoeffs(DFBout,parameters.deltaSize);

% Calculate scalar measures from cepstral coeffs:
output(1).ccMeasures = calculateHFCCmeasures(output(1).HFCC.cc);
output(1).DccMeasures = calculateHFCCmeasures(output(1).HFCC.Dcc);
output(1).DDccMeasures = calculateHFCCmeasures(output(1).HFCC.DDcc);

% Calculate scalar measures from FB output:
output(1).FBMeasures = calculateHFCCmeasures(FBout);
output(1).DFBMeasures = calculateHFCCmeasures(DFBout);
output(1).DDFBMeasures = calculateHFCCmeasures(DDFBout);

return;

% Bye!
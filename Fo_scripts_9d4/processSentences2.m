function [m,parameters] = processSentences2(m,fs,parameters)
% This function processes a struct of utterances of sentences for HFCC variation measures.  The
% following steps are followed:
% 1) Resample x to standard sampling rate (fsHFCC) for the sake of cepstral analysis,
% 2) High-pass filter to remove low-frequency noise,
% 3) Trim endpoints (silence regions) at beginning/end of signal using energy envelope,
% 4) Calculate active speech level and activity factor (ITU-T P.56 Method B),
% 5) Calculate cepstral coefficients (HFCC or MFCC) and delta/delta-delta coefficients,
% 6) Calculate CC/delta/delta-delta standard deviation sums.
%
% Input:
%    m -- struct of input data
%       .x -- Nx1 real vector, time series of utterance
%       .fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters used to process x [default values]
%       .fsHFCC -- real scalar, Hz, resampled rate of x [16000]
%       .HPF_order -- integer scalar, order of high-pass filter for removing low-freq. noise [10]
%       .HPF_fc -- real scalar, Hz, cutoff frequency of high-pass filter
%       .trimThreshold -- real scalar, dB, energy threshold for endpoint trimming [25]
%       .trimTimeConstant -- real scalar, sec, exp. moving average time constant [20e-3]
%       .calcASL -- logical scalar, true = calculate active speech level, false = don't [true]
%    Note: see HFCC.m for HFCC parameters.  Set parameters.fb to [] to design filter bank and DCT
%          matrix.
% Output:
%    m -- struct of output variables
%       .lenXtrimmed -- real scalar, sec, trimmed length of x
%       .ASL -- struct of active speech level function output (see ASL.m)
%       .activeSpeechLevel -- real scalar, dB DFS, active speech level of x
%       .activityFactor -- real scalar, unitless, ratio of active duration and total duration of x
%       .HFCC -- struct of HFCC output (see HFCC.m)
%       .ccMeasures -- struct of measures of HFCC matrix
%          .stDevSum -- real scalar, unitless, HFCC standard deviation sum measure
%          .varSum -- real scalar, unitless, HFCC variance sum measure
%          .euclDistMean -- real scalar, unitless, average Euclidean distance between adjacent
%                           frames
%          .absDistMean -- real scalar, unitless, average absolute distance between adjacent frames
%       .DccMeasures -- same as ccMeasures, but for delta HFCC features
%       .DDccMeasures -- same as ccMeasures, but for delta-delta HFCC features
%       .FBMeasures, .DFBMeasures, .DDFBMeasures -- same as CC measures, calculated from filter bank
%          output
%    parameters -- same as input, but with added default fields

% Mark Skowronski, May 29, 2013

% Set defaults for sentence preparations:
defaults(1).fsHFCC = 16000; % Hz
defaults(1).HPF_order = 10;
defaults(1).HPF_fc = 70; % Hz
defaults(1).trimThreshold = 25; % dB
defaults(1).trimTimeConstant = 20e-3; % sec
defaults(1).calcASL = true; 

% Check inputs:
if nargin<2,
   parameters = defaults;
else % add defaults fields to parameters if absent
   parameters = checkParameters(parameters,defaults);
end;

% Design Butterworth HPF:
[zHPF,pHPF,kHPF] = butter(parameters.HPF_order,parameters.HPF_fc/(parameters.fsHFCC/2),'high'); % returns poles/zeros/gain

% Process each WAV file:
for p=1:length(m),
   % Resample to fsHFCC:
   [num,den] = rat(parameters.fsHFCC/m(p).fs);
   x = resample(m(p).x,num,den);
   fs = parameters.fsHFCC;

   % HPF x using 2nd-order filtering:
   x = xFilter(x,zHPF,pHPF,kHPF);

   % Trim endpoints of x:
   x = trimEndpoints(x,fs,parameters.trimThreshold,parameters.trimTimeConstant);
   m(p).lenXtrimmed = length(x)/fs; % sec

   if parameters.calcASL,
      % Get active speech level:
      m(p).ASL = ASL(x,fs); % use default parameters
      m(p).activeSpeechLevel = m(p).ASL.activeSpeechLevel;
      m(p).activityFactor = m(p).ASL.activityFactor;
   end;

   % Get HFCCs and filter bank output:
   [m(p).HFCC,parameters,~,FBout] = HFCC(x,fs,parameters); % Designs filterbank and DCT matrix on first call
   
   % Calculate delta and delta-delta for FBout:
   DFBout = getDeltaCoeffs(FBout,parameters.deltaSize);
   DDFBout = getDeltaCoeffs(DFBout,parameters.deltaSize);

   % Calculate scalar measures from cepstral coeffs:
   m(p).ccMeasures = calculateHFCCmeasures(m(p).HFCC.cc);
   m(p).DccMeasures = calculateHFCCmeasures(m(p).HFCC.Dcc);
   m(p).DDccMeasures = calculateHFCCmeasures(m(p).HFCC.DDcc);
   
   % Calculate scalar measures from FB output:
   m(p).FBMeasures = calculateHFCCmeasures(FBout);
   m(p).DFBMeasures = calculateHFCCmeasures(DFBout);
   m(p).DDFBMeasures = calculateHFCCmeasures(DDFBout);
end;

return;

% Bye!
function [output,parameters,x,X] = processSentence(x,fs,parameters)
% This function processes an utterance of a sentence for HFCC variation measures and SWIPE measures
% of pitch.  The following steps are followed:
% 1) Resample x to fsHFCC for the sake of consistency,
% 2) [Optional] High-pass filter to remove low-frequency noise,
% 3) [Optional] Trim endpoints (silence regions) at beginning/end of signal using energy envelope,
% 4) [Optional] Calculate active speech level and activity factor (ITU-T P.56 Method B),
% 5) Calculate RMS value, in units of amplitude and dB DFS,
% 6) Calculate HFCCs and delta coefficients,
% 7) Calculate filter bank output (before cepstral transformation),
% 8) Calculate HFCC standard deviation sum and other measures,
% 9) [Optional] Calculate pitch and pitch strength and statistics using Aud-SWIPE'.
%
% Input:
%    x -- Nx1 real vector, time series of utterance
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters used to process x [defaults]
%    All measures:
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
%    HFCC parameters:
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
%    Active speech level parameters:
%       .calcASL -- logical scalar, true,[false] for calculating active speech level values
%    SWIPE parameters:
%       .calcSWIPE -- logical scalar, [true],false for calculating pitch and pitch strength values
%          and statistics using Aud-SWIPE'
%       .plim -- real 1x2 vector, Hz, [min max] frequency to consider for pitch [65 450]
%       .dt -- real scalar, sec, frame offset between adjacent analysis frames. Note: dt is the
%          inverse of frame rate [0.01] 
%       .dlog2p -- real scalar, octave, spacing between pitch candidates [1/24]
%       .dERBs -- real scalar, ERBrate, spacing between spectral estimates in ERBrate frequency
%          [0.10]
%       .woverlap -- real scalar in range [0,1], window overlap [0.5]
%       .sTHR -- real scalar in range [0,1], pitch strength threshold, pitch is NaN for frames with
%          pitch strength below threshold [0.1]
% Output:
%    output -- struct of output variables
%       .fs -- real scalar, Hz, resampled sampling rate, same as fsHFCC
%       .lenXtrimmed -- real scalar, sec, trimmed length of x
%       .tStart -- real scalar, sec, start time of trimmed x
%       .ASL -- struct of active speech level function output, [] if .calcASL = false
%       .activeSpeechLevel -- real scalar, dB DFS, active speech level of x
%       .activityFactor -- real scalar, unitless, ratio of active duration and total duration of x
%       .RMS, .RMSdB -- real scalar, root-mean-squared value of x after trimming and filtering
%       .HFCC -- struct of HFCC output (see HFCC.m)
%       Note: ASL struct and values are [] if calcASL is false.
%       .ccMeasures -- struct of output measures from cepstral coefficient matrix, HFCC.cc
%          .stDevVector -- Kx1 real vector, st. dev. of each ROW of X, K cepstral coefficients
%          .stDevSum -- real scalar, standard deviation along each ROW of X, summed
%          .varSum -- real scalar, variance along each ROW of X, summed
%          .euclDistMean -- real scalar, average of Euclidean distance between adjacent COLUMNS of X
%          .absDistMean -- real scalar, average of L1 distance between adjacent COLUMNS of X
%       .DccMeasures, .DDccMeasures -- same as .ccMeasures, for Dcc and DDcc matrices, respectively
%       .FBMeasures, .DFBMeasures, .DDFBMeasures, same as .ccMeasures, for filter bank matrix and
%          delta, delta-delta matrices, respectively
%       .ccStdSum -- same as .ccMeasures.stDevSum
%       .DccStdSum -- same as .DccMeasures.stDevSum
%       .DDccStdSum -- same as .DDccMeasures.stDevSum
%       .SWIPE -- struct of pitch, pitch strength, and stats ([] if calcSWIPE is false)
%          .P0time -- Gx1 real vector, sec, time points that begin analysis window for P0 and
%             P0strength
%          .P0 -- Gx1 real vector, Hz, pitch at each of G time points, P0(g) = NaN if
%             P0strength(g) < sTHR
%          .P0strength -- Gx1 real vector in range [0 1], pitch strength at each of G time points
%          .statsHz -- struct of stats calculated from P0 in Hz (NaN values ignored)
%             .mean, .median, .std, .min, .max -- real scalars, stats on non-NaN values in P0
%             .N -- integer scalar, number of non-NaN values in P0
%          .statsST -- same as .statsHz, for P0 in semitones
%          .statsPS -- same as .statsHz, for ALL pitch strength values (not just those above sTHR)
%    parameters -- same as input, but with added default fields
%    x -- Mx1 real vector, resampled and endpoint-trimmed version of input x, M samples
%    X -- spectrogram of x

% Mark Skowronski, September 18, 2013

% Set defaults for sentence preparations:
defaults(1).fsHFCC = 16000;
defaults(1).HPF_order = 10;
defaults(1).HPF_fc = 70;
defaults(1).trimThreshold = 25; % dB
defaults(1).trimTimeConstant = 0.02; % sec
defaults(1).ccMethod = 'HFCC';
defaults(1).windowType = 'Hamming';
defaults(1).windowSize = 0.0200; % sec
defaults(1).frameRate = 100;
defaults(1).FFTsize = 1024;
defaults(1).numFilters = 40;
defaults(1).numCoeffs = 13;
defaults(1).freqRange = [70 7500];
defaults(1).Efactor = 1;
defaults(1).deltaSize = 1;
defaults(1).applyPreemphasis = true;
defaults(1).preemphasisAlpha = 0.95;
defaults(1).applyCMS = true;
defaults(1).calcASL = false;
defaults(1).calcSWIPE = true;
defaults(1).plim = [65 450];
defaults(1).dt = 0.01; % sec
defaults(1).dlog2p = 0.0417; % octave
defaults(1).dERBs = 0.1000; % ERBrate
defaults(1).woverlap = 0.5;
defaults(1).sTHR = 0.1;

% Init output:
output(1).fs = defaults.fsHFCC;
output(1).lenXtrimmed = 0;
output(1).tStart = 0;
output(1).ASL = [];
output(1).activeSpeechLevel = [];
output(1).activityFactor = [];
output(1).RMS = [];
output(1).RMSdB = [];
output(1).HFCC = [];
output(1).ccMeasures = [];
output(1).DccMeasures = [];
output(1).DDccMeasures = [];
output(1).FBMeasures = [];
output(1).DFBMeasures = [];
output(1).DDFBMeasures = [];
output(1).ccStdSum = [];
output(1).DccStdSum = [];
output(1).DDccStdSum = [];
output(1).SWIPE = [];
X = [];

% Check inputs:
if nargin<3,
   parameters = defaults;
else % add defaults fields to parameters if absent
   parameters = checkParameters(parameters,defaults);
end;
if nargin<1 || isempty(x)
   return; % return empty output variables and default parameters
end;

% Resample to fsHFCC if necessary:
if parameters.fsHFCC~=fs,
   [num,den] = rat(parameters.fsHFCC/fs);
   x = resample(x,num,den);
end;
output(1).fs = parameters.fsHFCC; % Hz
fs = parameters.fsHFCC;

% HPF x using 2nd-order filtering if necessary:
if parameters.HPF_order>0 && parameters.HPF_fc>0,
   % Design and apply high-pass filter:
   [zHPF,pHPF,kHPF] = butter(parameters.HPF_order,parameters.HPF_fc/(fs/2),'high'); % returns poles/zeros/gain
   x = xFilter(x,zHPF,pHPF,kHPF);
end;

% Trim endpoints of x:
[x,t] = trimEndpoints(x,fs,parameters.trimThreshold,parameters.trimTimeConstant); % returns input x if time constant = 0
output(1).lenXtrimmed = length(x)/fs; % sec
output(1).tStart = t(1); % sec, start of x AFTER trimming

% Get active speech level if necessary:
if parameters.calcASL,
   output(1).ASL = ASL(x,fs); % use default parameters
   output(1).activeSpeechLevel = output(1).ASL.activeSpeechLevel;
   output(1).activityFactor = output(1).ASL.activityFactor;
else % create empty ASL outputs
   output(1).ASL = [];
   output(1).activeSpeechLevel = [];
   output(1).activityFactor = [];
end;

% Get RMS of (trimmed) x:
output(1).RMS = sqrt(mean(x.^2));
output(1).RMSdB = 20*log10(output.RMS);

% Get HFCCs:
[output(1).HFCC,parameters,X,FBout] = HFCC(x,fs,parameters); % Designs filterbank and DCT matrix on first call

% Calculate scalar measures from cepstral coeffs:
output(1).ccMeasures = calculateHFCCmeasures(output(1).HFCC.cc);
output(1).DccMeasures = calculateHFCCmeasures(output(1).HFCC.Dcc);
output(1).DDccMeasures = calculateHFCCmeasures(output(1).HFCC.DDcc);

% Calculate scalar measures from FB output:
DFBout = getDeltaCoeffs(FBout,parameters.deltaSize);
DDFBout = getDeltaCoeffs(DFBout,parameters.deltaSize);
output(1).FBMeasures = calculateHFCCmeasures(FBout);
output(1).DFBMeasures = calculateHFCCmeasures(DFBout);
output(1).DDFBMeasures = calculateHFCCmeasures(DDFBout);

% Legacy output fields:
output(1).ccStdSum = output(1).ccMeasures(1).stDevSum;
output(1).DccStdSum = output(1).DccMeasures(1).stDevSum;
output(1).DDccStdSum = output(1).DDccMeasures(1).stDevSum;

% Calculate pitch and pitch strength and statistics if necessary:
if parameters.calcSWIPE,
   % Calculate pitch and pitch strength:
   [P0,P0time,P0strength] = audswipep(x,fs,parameters.plim,parameters.dt,parameters.dlog2p,...
      parameters.dERBs,parameters.woverlap,parameters.sTHR);
   
   % Store:
   output(1).SWIPE(1).P0time = P0time;
   output(1).SWIPE(1).P0 = P0;
   output(1).SWIPE(1).P0strength = P0strength;
   
   % Calculate stats:
   output(1).SWIPE(1).statsHz = calculatePitchStats(P0); % Hz
   output(1).SWIPE(1).statsST = calculatePitchStats(f2ST(P0)); % semitones
   output(1).SWIPE(1).statsPS = calculatePitchStats(P0strength); % pitch strength
else % create empty SWIPE struct
   output(1).SWIPE = [];
end;

return;

function output = calculatePitchStats(f)
% This function calculates statistics for pitch frequencies f. NaN and inf values in f are ignored.

% Remove NaN and inf values in f:
f = f(~isnan(f) & ~isinf(f));

% Calculate statistics:
output(1).mean = mean(f);
output(1).median = median(f);
output(1).std = std(f);
output(1).min = min(f);
output(1).max = max(f);
output(1).N = length(f);

return;

% Bye!
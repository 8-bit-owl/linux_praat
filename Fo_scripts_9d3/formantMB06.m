function output = formantMB06(x,fs,parameters)
% This function is a wrapper for the formant estimation function formant_tracker_backend.m.  The
% input parameters and output variables are encapsulated in structs for convenience. The formant
% tracker requires 8 kHz sampling rate, so x is resampled to 8 kHz and returned in the output along
% with a vector of resampled time.
% Input:
%    x -- Nx1 real vector, time series data
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of control parameters [defaults]
%       .windowType -- string, ['Hamming'],'Hanning','Blackman','Rectangular' for spectral estimates 
%       .windowSize -- real scalar, sec, analysis window length [0.02]
%       .cutoffFreq -- real scalar, Hz, low-pass/high-pass cutoff frequency for log spectral density
%          ratio calculation, in range [700,1120] [700]
%       .logSpectralDensityRatioThreshold1 -- real scalar, unitless, unvoiced-to-voiced threshold,
%          in range [0.1,0.2] [0.2]
%       .logSpectralDensityRatioThreshold2 -- real scalar, unitless, voiced-to-unvoiced threshold,
%          in range [0.2,0.3] [0.3]
%       .autocorrThreshold -- real scalar, unitless, autocorrelation threshold for voicing
%          detection, in range [0.25,0.6] [0.4]
%       .RMSratio -- 1x4 real vector, dB, RMS thresholds for [F1,F2,F3,F4] adaptation
%          [-35,-40,-45,-50]
% Note: all input parameters except window type and size are adapted for each sample of input during
% estimation of formant frequencies.  All adapted variables except RMSratio vary in the specified
% ranges for each variable.
%
% Output:
%    output -- struct of output variables
%       .x -- Mx1 real vector, input resampled at 8 kHz, M samples
%       .time -- Mx1 real vector, sec, time points for samples in x
%       .fs -- real scalar, Hz, resampled frequency of 8 kHz
%       .F0 -- Mx1 real vector, Hz, fundamental frequency estimate
%       .F1 - .F4 -- Mx1 real vector, Hz, formant frequencies for F1 - F4
%       .voiceDecision -- Mx1 real vector, 0==unvoiced, 1==voiced
%       .gender -- Mx1 real vector, 0==male, 1==female, -1==previous

% Mark Skowronski, December 19, 2013

% Set parameter defaults:
defaults = promptFormantMB06Parameters([],true);

% Check inputs:
if nargin<3
   parameters = defaults;
else
   parameters = checkParameters(parameters,defaults);
end;
if nargin<2
   fs = 1;
end;
if nargin<1
   x = [];
end;

if ~isempty(x),
   % Resample x to 8 kHz:
   [aa,bb] = rat(8000/fs);
   x = resample(x,aa,bb);
   fs = 8000; % Hz
   time = [0:length(x)-1]'/fs; % sec, COLUMN vector

   % Translate input parameters (legacy variable names):
   lpc_window_size = round(fs*parameters.windowSize); % samples
   switch lower(parameters.windowType),
      case 'hamming'
         window = hamming(lpc_window_size);
      case 'hanning'
         window = hanning(lpc_window_size);
      case 'blackman'
         window = blackman(lpc_window_size);
      case 'rectangular'
         window = ones(lpc_window_size);
   end;
   window = repmat(window,1,4); % replicate for 4 formant channels
   Filter_Cutoff_init = parameters.cutoffFreq;
   Log_ratio_threshold1_init = parameters.logSpectralDensityRatioThreshold1;
   Log_ratio_threshold2_init = parameters.logSpectralDensityRatioThreshold2;
   Autocorrelation_Threshold_Level_init = parameters.autocorrThreshold;
   RMS_Ratio_F1 = parameters.RMSratio(1);
   RMS_Ratio_F2 = parameters.RMSratio(2);
   RMS_Ratio_F3 = parameters.RMSratio(3);
   RMS_Ratio_F4 = parameters.RMSratio(4);

   % Call formant tracker backend:
   [F,Voice,Pitch,Gender] = formant_tracker_backend(x,fs,lpc_window_size,window,Filter_Cutoff_init,...
      Log_ratio_threshold1_init,Log_ratio_threshold2_init,Autocorrelation_Threshold_Level_init,...
      RMS_Ratio_F1,RMS_Ratio_F2,RMS_Ratio_F3,RMS_Ratio_F4);
else
   % Init defaults:
   time = [];
   F = zeros(4,0); % empty vector with 4 rows
   Voice = [];
   Pitch = [];
   Gender = [];
end;

% Store output:
output(1).x = x;
output(1).time = time;
output(1).fs = fs;
output(1).F0 = Pitch;
output(1).F1 = F(1,:)'; % Hz, COLUMN vector
output(1).F2 = F(2,:)'; % Hz, COLUMN vector
output(1).F3 = F(3,:)'; % Hz, COLUMN vector
output(1).F4 = F(4,:)'; % Hz, COLUMN vector
output(1).voiceDecision = Voice'; % COLUMN vector
output(1).gender = Gender'; % 0=male, 1=female, -1=previous, COLUMN vector

return;
 
% Bye!
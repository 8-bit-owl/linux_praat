function parameters = promptFormantMB06Parameters(parameters,defaultsOnly)
% This function prompts for formantMB06 parameters or returns defaults according to the input flag.
% Input:
%    parameters -- scalar struct or empty, any input fields
%    defaultsOnly -- logical scalar, [false]: prompt user for parameter values, true: skip prompt
%       and return default values
% Output:
%    parameters -- scalar struct, same as input with added fields
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

% Set defaults:
defaults(1).windowType = 'Hamming';
defaults(1).windowSize = 20e-3; % sec
defaults(1).cutoffFreq = 700; % Hz
defaults(1).logSpectralDensityRatioThreshold1 = 0.2;
defaults(1).logSpectralDensityRatioThreshold2 = 0.3;
defaults(1).autocorrThreshold = 0.4;
defaults(1).RMSratio = [-35,-40,-45,-50]; % dB

% Check inputs:
if nargin<2,
   defaultsOnly = false; % prompt user for inputs
end;
if nargin<1 || isempty(parameters), % set defaults
   parameters = defaults;
else % Ensure all default fields exist in parameters:
   parameters = checkParameters(parameters,defaults);
end;

% Set parameter prompts and default values:
parameterBlock = {...
   'Window type [Hamming, Hanning, Blackman, Rectangular]',parameters.windowType;...
   'Window size, ms',num2str(parameters.windowSize*1e3);...
   'Cutoff frequency, Hz, range: [700,1120]',num2str(parameters.cutoffFreq);...
   'Log spectral density ratio threshold 1, range: [0.1,0.2]',num2str(parameters.logSpectralDensityRatioThreshold1);...
   'Log spectral density ratio threshold 2, range: [0.2,0.3]',num2str(parameters.logSpectralDensityRatioThreshold2);...
   'Auto-correlation threshold, range: [0.25,0.6]',num2str(parameters.autocorrThreshold);...
   'RMS threshold, dB, for F1',num2str(parameters.RMSratio(1));...
   'RMS threshold, dB, for F2',num2str(parameters.RMSratio(2));...
   'RMS threshold, dB, for F3',num2str(parameters.RMSratio(3));...
   'RMS threshold, dB, for F4',num2str(parameters.RMSratio(4))};

% Split into smaller groups:
displayString1 = parameterBlock(:,1);
valueString1 = parameterBlock(:,2);

% Prompt for parameter settings if necessary:
if ~defaultsOnly, % prompt user for values
   a = inputdlg(displayString1,'FORMANT PARAMETERS',[1 80],valueString1); % [] if cancel button
else % use default values without prompting
   a = []; % same as cancel button
end;

% Handle cancel button:
if isempty(a), % cancel button
   a = valueString1; % use initial values
end;

% Store input parameters:
parameters(1).windowType = a{1};
parameters(1).windowSize = str2num(a{2})/1000; % sec
parameters(1).cutoffFreq = str2num(a{3});
parameters(1).logSpectralDensityRatioThreshold1 = str2num(a{4});
parameters(1).logSpectralDensityRatioThreshold2 = str2num(a{5});
parameters(1).autocorrThreshold = str2num(a{6});
parameters(1).RMSratio = [str2num(a{7}),str2num(a{8}),str2num(a{9}),str2num(a{10})];

return;

% Bye!
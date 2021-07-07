function signal2moment(x,fs,parameters)
% This function calculates moments of a spectral estimate of x with sampling rate fs.
% The function performs the following steps:
% 1) Convert x to spectrum:
%    -- power spectral density, linear frequency scale
%    -- power spectral density, bark frequency scale
%    -- power spectral density, ERBrate frequency scale
%    -- excitation pattern, ERBrate frequency scale
%    -- specific loudness pattern, ERBrate frequency scale
% 2) Calculate moments from spectrum:
%    -- first moment (mean), second moment (variance), 3rd and 4th moments
%    -- skewness and kurtosis
% Input:
%    x -- Lx1 real vector, time domain signal to analyze, L samples
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of relevant parameters
%       .spectrumType -- ['PSDlinear'],'PSDbark','PSDERBrate','Excitation','SpecificLoudness'
%       .transducer -- string, name of acoustic transducer, ['freefield'],'er2', see
%                      outerEarConversion.m, applied to excitation and specific loudness only
%       .FFTsize -- integer scalar, length of FFT used to estimate spectrum of each frame [1024]
%       .windowType -- string of FFT window name: ['hamming'],'hanning','blackman','rectangular'
%       .windowSize -- real scalar, sec, duration of analysis window [20e-3]
%       .frameRate -- real scalar, frames/sec, number of analysis windows per second [100]
%       .fcRange -- 1x2 real vector, Hz, [min max] frequencies of excitation filter bank center
%                   frequencies [50 7500]
%       .ERBspacing -- real scalar, ERBrate units, spacing between adjacent center frequencies
%                      (equal spacing in ERBrate frequency space) [0.1]
%       .auditoryFilterType -- string, ['GCC'],'roex', name of auditory filter method to use to
%                              calculate excitation pattern

% Default parameters:
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

% Check inputs:
if nargin<3
   parameters = defaults;
end;

% Convert x to spectrum:
[X,f] = signal2spectrum(x,fs,parameters); % spectrum in each COLUMN

% Calculate moments from spectrum:
output = spectrum2moment(X,f);

return;

% Bye!
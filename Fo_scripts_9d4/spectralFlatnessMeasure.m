function [xSFM,t,f] = spectralFlatnessMeasure(x,fs,parameters)
% This function calculates the spectral flatness measure (SFM) for each frame of x.
% Inputs [defaults]:
%    x -- Nx1 real vector, time domain signal
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of SFM parameters
%       .FFTsize -- integer scalar, length of FFT used to estimate spectrum of each frame [1024]
%       .windowType -- string of FFT window name: ['hamming'],'hanning','blackman','rectangular'
%       .windowSize -- real scalar, sec, duration of analysis window [20e-3]
%       .frameRate -- real scalar, frames/sec, number of analysis windows per second [100]
% Outputs:
%    xSFM -- 1xT real vector, spectral flatness measure for each of T frames of x
%    t -- 1xT real vector, sec, T analysis frames, time of BEGINNING of each analysis frame
%    f -- Mx1 real vector, Hz, FFT bin frequencies

% Mark Skowronski, November 5, 2013

% Convert x to spectrogram (PSD):
[X,f,t] = signal2spectrum(x,fs,parameters);

% Skip DC bin:
X = X(2:end,:);

% Calculate geometric mean, arithmetic mean:
XgeoMean = exp(mean(log(X),1)); % mean along each COLUMN
XarithMean = mean(X,1);

% Calculate spectral flatness measure:
xSFM = XgeoMean./XarithMean;

return;

% Bye!
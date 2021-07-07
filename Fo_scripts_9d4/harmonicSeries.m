function x = harmonicSeries(fs,F0,T,slope,N)
% This function produces a harmonic series of a given spectral slope and fundamental frequency.
% Input:
%  fs -- real scalar, Hz, sampling rate of output
%  F0 -- real scalar, Hz, fundamental frequency of harmonic series
%  T -- real scalar, sec, duration of output
%  slope -- real scalar, dB/octave, rate of change of spectral peaks
%  N -- integer vector, harmonic indices to include in harmonic series, values in N >= 1, harmonics
%       above fs/2 ignored
% Output:
%  x -- Mx1 real vector, harmonic series output, signal scaled such that sin(F0) term has unity
%       amplitude, M = round(fs*T)

% Mark Skowronski, November 30, 2011

% Init output:
M = round(fs*T); % samples
x = zeros(M,1);

% Check N:
maxN = ceil((fs/2)/F0)-1; % maximum harmonic index < fs/2
N = N(N<=maxN & N>0); % ignore harmonic indices above fs/2 and <= 0

% Create time vector:
t = [0:M-1]'/fs; % sec, COLUMN vector

% Add each harmonic to x:
for p=1:length(N),
   % Determine frequency of harmonic component:
   F0p = N(p)*F0; % Hz
   
   % Determine amplitude of harmonic component:
   numOctaves = log(F0p/F0)/log(2); % number of octaves from F0 to F0p
   AdB = numOctaves*slope; % dB, amplitude of harmonic component, assuming F0 harmonic at 0 dB
   A = 10^(AdB/20); % amplitude units
   
   % Add harmonic component to x:
   x = x + A*sin(2*pi*F0p*t);
end;

return;

% Bye!
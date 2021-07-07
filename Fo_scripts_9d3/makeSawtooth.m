function x = makeSawtooth(T,F0,fs,spectralSlopeFlag,K,phaseFlag)
% This function creates a sawtooth wave with harmonics up to the Nyquist rate.  The sawtooth is
% constructed with scaled sinusoids up to the Nyquist frequency (unaliased).
% Input:
%  T -- real scalar, sec, duration of signal
%  F0 -- real scalar, Hz, fundamental frequency of sawtooth
%  fs -- real scalar, Hz, sampling rate of signal
%  spectralSlopeFlag -- integer scalar, indicates signal type
%     1: sawtooth wave -- harmonics attenuate 6 dB/octave (default)
%     2: pulse train -- harmonics with equal intensity
%  K -- integer scalar, # harmonics to include in harmonic complex
%  phaseFlag -- integer scalar, indicates relative phases of harmonics
%     1: 0 phase for all harmonics (default)
%     2: uniform random phase in range [0,2*pi] radians
% Output:
%  x -- Nx1 real vector, sawtooth signal, N=round(T*fs) samples

% Mark Skowronski, November 8, 2011

Kmax = floor((fs/2)/F0); % largest harmonic <= Nyquist rate

% Check inputs:
if nargin<6,
   phaseFlag = 1;
end;
if nargin<5,
   K = Kmax;
else
   K = min(K,Kmax);
end;
if nargin<4,
   spectralSlopeFlag = 1;
end;

% Construct x as harmonic series:
N = round(T*fs); % samples, duration of x
x = zeros(N,1); % init
t = [0:N-1]'/fs; % sec, COLUMN vector
for p=1:K, % for each harmonic
   switch spectralSlopeFlag
      case 1 % sawtooth
         A = -2/pi/p;
      case 2 % pulse train
         A = 1/K;
   end;
   switch phaseFlag
      case 1 % zero phase
         ph = 0;
      case 2 % uniform random
         ph = rand*2*pi;
   end;
   x = x + A*sin(2*pi*p*F0*t+ph);
end;


return;


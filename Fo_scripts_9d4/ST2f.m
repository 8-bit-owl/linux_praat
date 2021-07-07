function f = ST2f(ST)
% This function converts frequency in semitones to Hz.
% Input:
%    ST -- Nx1 real vector, semitones, frequencies to convert to Hz
% Output:
%    f -- Nx1 real vector, Hz
%
% Reference: R. W. Young (1939). "Terminology for logarithmic frequency units," JASA, 11(1), pp.
% 134-139

% Mark Skowronski, September 18, 2013

% Set reference:
A4 = 440; % Hz
C5 = A4*2^(3/12); % Hz, C above A440
C0 = C5*2^(-5); % Hz, C0 ~ 16.352 Hz = 0 semitones

% Calculate semitones:
f = exp(ST/12*log(2))*C0;

return;

% Bye!
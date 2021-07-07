function ST = f2ST(f)
% This function converts frequency in Hz to semitones.
% Input:
%    f -- Nx1 real vector, Hz, frequencies to convert to semitones
% Output:
%    ST -- Nx1 real vector, semitones
%
% Reference: R. W. Young (1939). "Terminology for logarithmic frequency units," JASA, 11(1), pp.
% 134-139

% Mark Skowronski, September 18, 2013

% Set reference:
A4 = 440; % Hz
C5 = A4*2^(3/12); % Hz, C above A440
C0 = C5*2^(-5); % Hz, C0 ~ 16.352 Hz = 0 semitones

% Calculate semitones:
ST = 12*log(f/C0)/log(2);

return;

% Bye!
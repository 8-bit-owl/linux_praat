function fERBrate = f2ERBrate(f)
% This function converts frequency to ERB-rate according to Glasberg and Moore (1990).
% Input:
%    f -- Nx1 real vector, Hz, frequencies to convert to ERB-rate
% Output:
%    fERBrate -- Nx1 real vector, Hz/Hz, ERB-rate of f
%
% Reference: Glasberg, B. R. and B. C. J. Moore, "Derivation of auditory filter shapes from
% notched-noise data," Hearing Research, vol. 47, pp. 103:138, 1990
%

% Mark Skowronski, December 4, 2012

% Set scalars (p. 132, FORTRAN code):
c1 = 24.673;
c2 = 4.368;

% Calculate fERB:
fERBrate = 1000/(c1*c2)*log(c2/1000*f+1); % Hz/Hz, Equation 4, p. 114

return;

% Bye!
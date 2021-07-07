function f = ERBrate2f(fERBrate)
% This function converts ERB-rate to linear frequency according to Glasberg and Moore (1990).
% Input:
%    fERBrate -- Nx1 real vector, Hz/Hz, ERB-rate of N frequencies
% Output:
%    f -- Nx1 real vector, Hz, linear frequency of ERB-rate
%
% Reference: Glasberg, B. R. and B. C. J. Moore, "Derivation of auditory filter shapes from
% notched-noise data," Hearing Research, vol. 47, pp. 103:138, 1990
%

% Mark Skowronski, December 4, 2012

% Set scalars (p. 132, FORTRAN code):
c1 = 24.673;
c2 = 4.368;

% Calculate fERB:
f = (exp((c1*c2)/1000*fERBrate)-1)*1000/c2; % Hz, inverse of Eq. 4, p. 114

return;

% Bye!
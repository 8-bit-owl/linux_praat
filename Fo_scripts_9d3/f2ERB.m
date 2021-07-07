function fERB = f2ERB(f)
% This function converts frequency to ERB frequency according to Glasberg and Moore (1990).
% Input:
%    f -- Nx1 real vector, Hz, frequencies to convert to ERB
% Output:
%    fERB -- Nx1 real vector, Hz, ERB frequencies of f
%
% Reference: Glasberg, B. R. and B. C. J. Moore, "Derivation of auditory filter shapes from
% notched-noise data," Hearing Research, vol. 47, pp. 103:138, 1990
%
% Note: The linear function converting f to fERB replaces the quadratic function from Moore and
% Glasberg (1983).

% Mark Skowronski, December 3, 2012

% Set scalars (p. 132, FORTRAN code):
c1 = 24.673;
c2 = 4.368;

% Calculate fERB:
fERB = c1*(c2*f/1000+1); % Hz, Equation 3, p. 114

return;

% Bye!
function f = bark2f(z)
% This function converts frequencies z in Bark to Hz.
% Input:
%    z -- Mx1 real vector, Bark, frequencies
% Output:
%    f -- Mx1 real vector, Hz, converted frequencies
%
% Note: No closed-form expression exists converting z to f, so reference vectors are constructed for
% interpolation. The reference frequency vector includes integer values between 0 and 50 kHz.
%
% Reference: Zwicker and Terhardt, JASA 68(5), pp. 1523-1525, Nov. 1980

% Mark Skowronski, May 20, 2013

% Construct references for interpolation:
fRef = [0:50e3]; % Hz, 1-Hz spacing
zRef = f2bark(fRef); % Bark

% Interpolate z using reference to estimate f (no closed-form expression):
f = interp1(zRef,fRef,z,'linear','extrap');

return;

% Bye!
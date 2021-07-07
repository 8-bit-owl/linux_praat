function [W,f] = cosPowerFunction(f,fr,A,BW)
% This function creates a cosine-power auditory function with resonant frequency fr, asymptote level
% A, and bandwidth BW.  The returned function is in power units.
% W(f) = (1-Ap)/2*(1-cos(pi*(f/fr)^pl))+Ap, f<=fr
%      = 1/2*(1-cos(pi*(f/fr)^pu)), f>fr
%    Ap -- A, in power units
%    pl,pu -- bandwidth control factor for lower, upper skirts, respectively
% Input [defaults]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate function [0:10000]
%    fr -- real scalar, Hz, resonant frequency of function [1000]
%    A -- real scalar, dB, asymptote floor of function, W(f=0) = A dB [-40]
%    BW -- real scalar, Hz, half-power bandwidth of function [100]
% Output:
%    W -- Nx1 real vector, squared magnitude spectrum of auditory function
%    f -- Nx1 real vector, Hz, same as input; otherwise, default f values.

% Mark Skowronski, March 25, 2013

% Check inputs:
if nargin<4,
   BW = 100; % Hz
end;
if nargin<3,
   A = -40;
end;
if nargin<2
   fr = 1000; % Hz
end;
if nargin<1
   f = [0:10000]'; % Hz, COLUMN vector
end;

% Ensure f is COLUMN vector:
f = f(:);

% Convert A to power units:
Ap = 10^(A/10); % power units

% Compute lower and upper exponent terms:
[pl,pu] = getExponent(fr,Ap,BW);

% Compute lower and upper warped frequencies:
[flhat,fuhat] = getFhat(f,fr,pl,pu);

% Compute lower and upper W:
Wl = (1-Ap)/2*(1-cos(pi*flhat))+Ap;
Wu = 1/2*(1-cos(pi*fuhat)); % no lower asymptote

% Concatenate lower and upper W:
W = [Wl;Wu];

return;

function [flhat,fuhat] = getFhat(f,fr,pl,pu)
% This function computes warped frequency fhat below and above fr.
% fhat(0) = 0
% fhat(f) monotonically increases from 0 to 1 for f in [0,fr]
% fhat(fr) = 1
% fhat(f) monotonically decreases from 1 to 0 for f in [fr,inf)
% fhat(f-->inf) --> 0

% Separate f below and above fr:
fl = f(f<=fr);
fu = f(f>fr);

% Calculate fhat below and above fr:
flhat = (fl/fr).^pl;
fuhat = (fu/fr).^pu;

return;

function [pl,pu] = getExponent(fr,Ap,BW)
% This function calculates the lower and upper exponent terms for fhat.

% Calculate pre-terms:
A1 = acos(Ap/(1-Ap))/pi; % ~1/2 for A<<1
fcl = max(1,fr-BW/2); % lower half-power cutoff frequency, minimum: 1 Hz
fcu = fr+BW/2; % upper half-power cutoff frequency

% Calculate pl and pu:
pl = log(A1)/(log(fcl/fr));
pu = log(1/2)/(log(fcu/fr)); % no A term for upper skirt, so use log(1/2)

return;

% Bye!

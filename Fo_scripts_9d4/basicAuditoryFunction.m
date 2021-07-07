function [W,f] = basicAuditoryFunction(f,fr,A,BW,fType)
% This function creates a basic auditory function with resonant frequency fr, asymptote level A, and
% bandwidth BW.  The returned function is in power units.
% W(f) = (1-Ap)/2*(1-cos(pi*fhat))+Ap
%    Ap -- A, in power units
%    fhat -- warped frequency, exponential or power
%        exp: fhat(f) = (exp(p*f/fr)-1)/(exp(p)-1)
%        power: fhat(f) = (f/fr)^p
%        p -- bandwidth control factor
% Input [defaults]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate function [0:10000]
%    fr -- real scalar, Hz, resonant frequency of function [1000]
%    A -- real scalar, dB, asymptote floor of function, W(f=0) = A dB [-40]
%    BW -- real scalar, Hz, half-power bandwidth of function [100]
%    fType -- string, ['power'],'exp', frequency-warping function
% Output:
%    W -- Nx1 real vector, squared magnitude spectrum of auditory function
%    f -- Nx1 real vector, Hz, same as input; otherwise, default f values.

% Mark Skowronski, March 25, 2013

% Check inputs:
if nargin<5,
   fType = 'power';
end;
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
[pl,pu] = getExponent(fr,Ap,BW,fType);

% Compute lower and upper warped frequencies:
[flhat,fuhat] = getFhat(f,fr,pl,pu,fType);

% Compute lower and upper W:
Wl = (1-Ap)/2*(1-cos(pi*flhat))+Ap;
Wu = 1/2*(1-cos(pi*fuhat)); % no lower asymptote

% Concatenate lower and upper W:
W = [Wl;Wu];

return;

function [flhat,fuhat] = getFhat(f,fr,pl,pu,fType)
% This function computes warped frequency fhat below and above fr.
% fhat(0) = 0
% fhat(f) monotonically increases from 0 to 1 for f in [0,fr]
% fhat(fr) = 1
% fhat(f) monotonically decreases from 1 to 0 for f in [fr,inf)
% fhat(f-->inf) --> 0

% Separate f below and above fr:
fl = f(f<=fr);
fu = f(f>fr);

% Calculate fhat below and above fr depending on fType:
switch lower(fType),
   case 'exp'
      flhat = (exp(pl*fl/fr)-1)/(exp(pl)-1); % lower frequency
      fuhat = exp(pu*(fu/fr-1)); % upper frequency
   otherwise % 'power'
      flhat = (fl/fr).^pl;
      fuhat = (fu/fr).^pu;
end;

return;

function [pl,pu] = getExponent(fr,Ap,BW,fType)
% This function calculates the lower and upper exponent terms for fhat depending on fType.

% Calculate pre-terms common to all fTypes:
A1 = acos(Ap/(1-Ap))/pi; % ~1/2 for A<<1
fcl = fr-BW/2; % lower half-power cutoff frequency
fcu = fr+BW/2; % upper half-power cutoff frequency

switch lower(fType),
   case 'exp'
      pl = log(A1)/(fcl/fr-1); % approximation, exp(p)>>1
      pu = log(A1)/(fcu/fr-1); % approximation, exp(p)>>1
   otherwise % 'power'
      pl = log(A1)/(log(fcl)-log(fr));
      pu = log(A1)/(log(fcu)-log(fr));
end;

return;

% Bye!

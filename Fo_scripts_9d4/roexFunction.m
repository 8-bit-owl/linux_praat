function [W,f] = roexFunction(f,fc,pu,pl)
% This function creates a rounded-exponential (roex) function used in the excitation model of
% Glasberg & Moore 1990.
% Roex function: h(g) = (1+p*abs(g))*exp(-p*abs(g)), g = (f-fc)/fc, p may be different for upper and
% lower tails of function on either side of fc.
% Input [defaults]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate roex function [0:2000]
%         Note: if f is empty, the function adds defaults to parameters and returns.
%    fc -- real scalar, Hz, center frequency of function [1000]
%    pu,pl -- real scalar, upper and lower p-values of roex [30,30]
% Output:
%    W -- Nx1 real vector, squared magnitude spectrum of auditory filter
%    f -- Nx1 real vector, Hz, same as input; otherwise, default f values.

% Mark Skowronski, January 18, 2013

% Check inputs:
if nargin<4,
   pl = 30;
end;
if nargin<3,
   pu = 30;
end;
if nargin<2
   fc = 1000; % Hz
end;
if nargin<1
   f = [0:2000]'; % Hz
end;

% Ensure f is COLUMN vector:
f = f(:);

% Split f into upper and lower tails:
fl = f(f<fc);
fu = f(f>=fc);

% Calculate roex for each tail:
Wl = getRoex(fl,fc,pl);
Wu = getRoex(fu,fc,pu);

% Combine:
W = [Wl;Wu];

return;


function W = getRoex(f,fc,p)
% This function calculates a roex function for single p value.

% Calculate formalized frequency:
g = abs((f-fc)/fc);

% Calculate roex function:
PG = p*g;
W = (1+PG).*exp(-PG);

return;

% Bye!
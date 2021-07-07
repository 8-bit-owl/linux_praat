function [x,t] = trimEndpoints(x,fs,trimThreshold,trimTimeConstant)
% This function trims the endpoints of x.  An exponential moving average of the power in x is
% calculated, and the beginning/end of x less than the threshold away from the peak averaged power
% are trimmed.
% Input [defaults]:
%    x -- Nx1 real vector, input signal to be trimmed, typically a hand-labeled utterance with
%       arbitrary silence at beginning/end of x
%    fs -- real scalar, Hz, sampling rate of x
%    trimThreshold -- real scalar, dB, beginning/end of x below threshold is trimmed [25]
%    trimTimeConstant -- real scalar, sec, time constant used to calculate energy contour [0.02]
%    Note: set time constant to 0 to skip endpoint trimming.
% Output:
%    x -- Mx1 real vector, M<=N
%    t -- Mx1 real vector, sec, time vector associated with output x, referenced to input x

% Check inputs:
if nargin<4
   trimTimeConstant = 0.02; % sec
end;
if nargin<3
   trimThreshold = 25; % dB
end;

% Create time vector:
t = [0:length(x)-1]'/fs; % sec

% Skip endpoint trimming if necessary:
if trimTimeConstant==0,
   return; % exit function without modifying x
end;

% Construct EMA filter coefficient:
g = exp(-1/(trimTimeConstant*fs));

% Filter power in x:
xPow = filter(1,[1 -g],x.^2); % arbitrary gain, since trimming is relative to max

% Convert xPow to dB:
xPow = 10*log10(xPow); % dB

% Find values in xPow within trimThreshold of max(xPow):
xPowMax = max(xPow);
xPowRange = find(xPow>=(xPowMax-trimThreshold));

% Trim x, t:
x = x(xPowRange(1):xPowRange(end));
t = t(xPowRange(1):xPowRange(end));

return;

% Bye!
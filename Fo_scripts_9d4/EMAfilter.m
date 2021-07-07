function y = EMAfilter(x,fs,T)
% This function performs exponential moving average (EMA) filtering on x.  For input x(n) and output
% y(n):
%    y(1) = x(1), initialization
%    y(n) = a*y(n-1) + (1-a)*x(n), n>1
% where a = exp(-1/(T*fs)) and T is the time constant of the filter.  As T goes to zero, then a goes
% to zero, and the filter output is simply the current input sample: no previous input samples
% affect the output.  As T goes to infinity, then a goes to one, and the filter output approaches a
% long-term average of the input: the current input sample has no effect on the (infinity) long-term
% average.
%
% Input:
%    x -- Nx1 real vector, time series signal to be filtered, length N samples
%    fs -- real scalar, Hz, sampling rate of x
%    T -- real scalar, sec, time constant of EMA filter (if T<=0, then y=x is returned)
% Output:
%    y -- Nx1 real vector, EMA filtered version of x

% Mark Skowronski, October 3, 2012

% Calculate filter coefficient from T:
if T<=0, % special case
   y = x;
   return;
end;
a = exp(-1/(T*fs));

% Init output:
y = zeros(size(x));
y(1) = x(1);

% Iterate output:
zi = y(1)*a; % initial condition for filter function
y(2:end) = filter([1-a],[1,-a],x(2:end),zi);
%for n=2:length(x),
%   y(n) = a*y(n-1) + (1-a)*x(n);
%end;

return;

% Bye!
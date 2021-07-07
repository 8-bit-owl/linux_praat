function [output,parameters] = ASL(x,fs,parameters)
% This function measures active speech level according to ITU P.56 Method B.
% Input
%  x -- Nx1 real vector, speech signal, abs(x) <= 1
%  fs -- integer scalar, Hz, sampling rate of x
%  parameters -- (optional) struct of ASL parameters
%     .f -- integer scalar, Hz, sampling rate of speech signal envelope (default: 694)
%     .T -- real scalar, sec, exponential moving average time constant (default: 0.03)
%     .H -- real scalar, sec, hangover time (default: 0.2)
%     .M -- real scalar, dB, margin (default: 15.9)
% Output
%  output -- struct of output values
%     .totalTime -- real scalar, sec, duration of decimated x
%     .longTermPower -- real scalar, volts^2, power of decimated x
%     .activePowerEstimate -- 1xJ real vector, volts^2, active power estimate of decimated x for
%        each of J thresholds with respect to reference voltage r
%     .L -- real scalar, dB, long-term level of decimated x
%     .A -- 1xJ real vector, dB, active-level estimate of decimated x for each of J thresholds with
%        respect to reference voltage r
%     .C -- 1xJ real vector, dB, thresholds used to calculate A
%     .qAll -- 1xn real vector, envelope of decimated x of length n
%     .c -- 1xJ real vector, thresholds used in comparison with qAll
%     .a -- 1xJ integer vector, active-level counts for each threshold
%     .n -- integer scalar, length of decimated x
%     .s -- real scalar, sum of samples of decimated x
%     .sq -- real scalar, sum of squared samples of decimated x
%     .activeSpeechLevel -- real scalar, dB, active speech level of input
%     .activeSpeechLevelThreshold  -- real scalar, dB, threshold used to determine ASL
%     .activityFactor -- real scalar, fraction of total time speech is active, in range [0,1]
%
% Constants:
%  r -- reference voltage = 1 volt
%  v -- D/A conversion factor = 1 volt/unit

% Mark Skowronski, March 15, 2010

% Check input:
if nargin<3,
   parameters = struct([]);
   parameters(1).f = 694; % Hz, sampling rate of envelope
   parameters(1).T = 0.03; % sec, exponential moving average time constant
   parameters(1).H = 0.2; % sec, hangover time
   parameters(1).M = 15.9; % dB, margin
end;

% Calculate terms from parameters:
t = 1/parameters.f; % sec, sampling interval of envelope
g = exp(-t/parameters.T);
I = ceil(parameters.H/t); % samples, hangover, rounded up

% Set D/A conversion constant and reference voltage:
v = 1; % volt/unit
r = 1; % volt, reference

% Create thresholds:
c = 10.^([-100:.1:0]/20); % -100 to 0 dB range, 0.1 dB steps
a = zeros(size(c)); % activity count
h = I*ones(size(c)); % hangover count

% Downsample x:
xd = x(round([1:fs/parameters.f:length(x)])); % disregard aliasing

% Process 1:
n = length(xd);
s = sum(xd);
sq = sum(xd.^2);

% Process 2 and activity/hangover count update:
qAll = zeros(1,n);
p = (1-g)*abs(xd(1));
q = (1-g)*p; % envelope
[a,h] = updateAH(a,h,q,c,I);
qAll(1) = q;
for i=2:n,
   % Update envelope:
   p = g*p + (1-g)*abs(xd(i));
   q = g*q + (1-g)*p;
   qAll(i) = q;
   
   % Update counts:
   [a,h] = updateAH(a,h,q,c,I);
end;

% Calculate output terms:
output = struct([]);
output(1).totalTime = n*t; % sec
output.longTermPower = sq/n*v^2; % volts^2
output.activePowerEstimate = sq./(a+1e-308)*v^2; % volts^2, add small offset to avoid division by zero
output.L = 10*log10(output.longTermPower) - 20*log10(r); % dB, long-term level
output.A = 10*log10(output.activePowerEstimate) - 20*log10(r); % dB, active-level estimate
output.C = 20*log10(c*v) - 20*log10(r); % dB, threshold
output.qAll = qAll; % envelope
output.c = c;
output.a = a;
output.n = n;
output.s = s;
output.sq = sq;

% Find active speech level:
[junk,minIndex] = min(abs((output.A - output.C) - parameters.M)); % index of threshold closest to A-C = M
output.activeSpeechLevel = output.A(minIndex);
output.activeSpeechLevelThreshold = output.C(minIndex);
output.activityFactor = a(minIndex)/n; % between 0 and 1

return;



function [a,h] = updateAH(a,h,q,c,I)
% This function updates the activity and hangover counts.

for j=1:length(c), % for each threshold
   if q>=c(j), % at or above threshold
      a(j) = a(j)+1;
      h(j) = 0;
   elseif h(j)<I, % below threshold but in hangover region
      a(j) = a(j)+1;
      h(j) = h(j)+1;
   end;
end;

return;

% Bye!
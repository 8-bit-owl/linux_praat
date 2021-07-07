function [t,F0,rms,zc,cp,cpp,parameters] = CPP1994(x,fs,parameters)
% This function calculates cepstral peak prominence (CPP) according to Hillenbrand et al. (1994).
% The function tries to duplicate the results from cpps.exe.
% Input:
%    x -- Lx1 real vector, time domain signal to analyze
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters
%       .windowLength -- integer scalar, samples, length of analysis window [512]
%       .shift -- real scalar, ms, frame offset of adjacent analysis frames [2]
%       .lower -- real scalar, Hz, lower limit to search for cepstral peak [60]
%       .upper -- real scalar, Hz, upper limit to search for cepstral peak [300]
%       .FFTsize -- integer scalar, length of FFT of each window [2^13]
% Output:
%    t -- Kx1 real vector, ms, time of analysis frame, K total frames
%    F0 -- Kx1 real vector, Hz, fundamental frequency according to cepstral peak lag
%    rms -- Kx1 real vector, dB, RMS level of frame 10 ms in duration (ind. of windowLength)
%    zc -- Kx1 integer vector, zero-crossing count of frame
%    cp -- Kx1 real vector, dB, cepstral peak before normalization
%    cpp -- Kx1 real vector, dB, cepstral peak after normalization
%    parameters -- same as input, or set to defaults if no input
%    
% Reference:
% Hillenbrand, Cleveland, and Erickson, "Acoustic Correlates of Breathy Vocal Quality," JSHR, vol.
% 37, pp. 769-778, Aug. 1994

% Mark Skowronski, June 25, 2013

defaults(1).windowLength = 512;
defaults(1).shift = 2;
defaults(1).lower = 60;
defaults(1).upper = 300;
defaults(1).FFTsize = 2^13;

% Check inputs:
if nargin<3
   parameters = defaults;
end;

% Process each frame:
p0 = 1; % sample, pointer to first sample in frame
k = 0; % frame index
T10 = round(fs*10e-3); % samples, duration of 10-ms window for RMS calculation
q1 = ceil(fs/1000+1); % samples, index into cepstrum C
while p0-1+max(parameters.windowLength,T10) <= length(x)
   % Get frame and 10-ms frame, update counters:
   k = k+1; % update frame index
   tOffset = round((k-1)*parameters.shift/1e3*fs); % samples, shift of current frame
   tIndex = [1:parameters.windowLength] + tOffset; % samples, index into x of current frame
   t10Index = [1:T10] + tOffset; % samples, index into x of 10-ms frame
   xFrame = x(tIndex);
   x10Frame = x(t10Index);
   t(k,1) = (p0-1)/fs*1e3; % ms
   p0 = 1 + round((k)*parameters.shift/1e3*fs); % samples, index into NEXT frame of x
   
   % Get RMS values from 10-ms frame:
   rms(k,1) = 10*log10(mean(x10Frame.^2)); % dB
   
   % Get zero-crossing count in frame:
   zc(k,1) = sum(xFrame(1:end-1)<=0 & xFrame(2:end)>0); % positive crossing count
   
   % Get real cepstrum of frame:
   X = 20*log10(abs(fft(xFrame))+1e-300); % dB, add small offset to avoid log(0)
%   X = adjustLogSpectrum(X,peakFactor);
   c = ifft(X,'symmetric'); % treat X as conjugate symmetric (c real)
   C = 20*log10(abs(c)); % dB
   
   % Determine limits over which to search for peak in C and to perform cepstral baseline regression:
   tRange = [ceil(fs/parameters.upper)+1,floor(fs/parameters.lower)+1]; % samples, index into C
   tRange = [max(1,tRange(1)),min(length(C)/2,tRange(2))]; % limit to integer indices into C
   
   % Get cepstral peak:
   [cp(k,1),CmaxIndex] = max(C(tRange(1):tRange(2)));
   CmaxIndex = CmaxIndex+tRange(1)-1; % index into C

   % Find regression of C in tRange:
   [m,b] = calcRegress(C(q1:end/2)); % [slope, y-intercept] over domain [1:N]
   
   % Calculate baseline at CmaxIndex:
   Cbaseline = m*(CmaxIndex-(q1-1))+b; % from regression
   cpp(k,1) = cp(k,1) - Cbaseline; % normalize with baseline value
   
   % Calculate F0:
   F0(k,1) = fs/(CmaxIndex-1); % Hz, convert index to lag with '-1'
end;

return;

function [yStar,xStar] = quadInterp(y)
% This function fits a quadratic spline to the 3 points in y(x), x=[-1,0,1], and returns the peak of
% the spline, yStar, at xStar.
% Spline: yHat = a0 + a1*x + a2*x^2

% Create mixing matrix for coefficients [a0,a1,a2]:
R = [1 -1 1;1 0 0;1 1 1];

% Solve for coefficients:
a = R\y(:); % equivalent: a = inv(R)*y(:), ensure y is COLUMN vector

% Solve for xStar:
xStar = -a(2)/(2*a(3));

% Solve for yStar:
yStar = a(1)+a(2)*xStar+a(3)*(xStar^2);

return;

function [m,b] = calcRegress(y)
% This function calculates a least-squares solution to the regression of y(x),
% x=[1:N], N=length(y).
% Input:
%    y -- Nx1 real vector, signal to regress
% Output:
%    m -- real scalar, slope of regression line
%    b -- real scalar, y-intercept of regression line

% Create x vector:
N = length(y);
x = [[1:N]',ones(N,1)]; % include COLUMN of 1s for y-intercept term

% Create inner and outer product terms:
R = x'*x; % inner product term
p = x'*y(:); % outer product term, ensure y is COLUMN vector

% Solve:
B = R\p; % equivalent: inv(R)*p
m = B(1); % slope
b = B(2); % y-intercept

return;

function X = adjustLogSpectrum(X,peakFactor)
% This function adjusts the log spectrum X to minimize sensitivity of IFFT to large negative log
% spectrum values

% Shift X:
X = X-max(X)+peakFactor; % X peak at peakFactor

% Truncate:
X(X<0) = 0;

return;

% Bye!
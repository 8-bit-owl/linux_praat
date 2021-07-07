function [P,F0,C] = CPP(x,fs,fRange,fftSize)
% This function calculates cepstral peak prominence (CPP) according to Hillenbrand et al. (1994).
% Input [defaults]:
%    x -- Lx1 real vector, time domain signal to analyze
%    fs -- real scalar, Hz, sampling rate of x
%    fRange -- 1x2 real vector, Hz, [min max] frequency to search for cepstral peak [50 500]
%    fftSize -- integer scalar, FFT size to use on x [2^15]
% Output:
%    P -- real scalar, dB, cepstral peak prominence value
%    F0 -- real scalar, Hz, frequency of cepstral peak in fRange
%    C -- Lx1 real vector, dB, real cepstrum of x
%    
% Reference:
% Hillenbrand, Cleveland, and Erickson, "Acoustic Correlates of Breathy Vocal Quality," JSHR, vol.
% 37, pp. 769-778, Aug. 1994

% Mark Skowronski, June 11, 2013

% Check inputs;
if nargin<4
   fftSize = 2^15;
end;
if nargin<3
   fRange = [50 500]; % Hz, range of F0 to search for cepstral peak
end;

% Spectrum magnitude:
Xabs = abs(fft(x(:),fftSize)); % COLUMN vector

% Smooth spectrum:
Hsmooth = [1 2 1]'/2; % COLUMN vector
Xabs = conv(Xabs,Hsmooth,'same');

% Log spectrum:
X = log(Xabs);

% Zero mean:
X = X-mean(X);

% Real cepstrum:
c = ifft(X,'symmetric'); % treat X as conjugate symmetric (c real)
C = 20*log10(abs(c)); % dB

% Determine limits over which to search for peak in C and to perform cepstral baseline regression:
tRange = [ceil(fs/fRange(2))+1,floor(fs/fRange(1))+1]; % samples, index into C
tRange = [max(1,tRange(1)),min(fftSize/2,tRange(2))]; % limit to integer indices into C

% Find peak of C in tRange:
CRange = C(tRange(1):tRange(2));
[Cmax,CmaxIndex] = max(CRange); % index into CRange in domain [1:N], N=length(CRange)

% Adjust Cmax, CmaxIndex using quadratic interpolation:
%[Cmax,deltaIndex] = quadInterp(C([CmaxIndex-1:CmaxIndex+1]+tRange(1)-1)); % index into C in case CmaxIndex=1
%CmaxIndex = CmaxIndex+deltaIndex;

% Find regression of C in tRange:
[m,b] = calcRegress(CRange); % [slope, y-intercept] over domain [1:N]

% Calculate cepstral peak prominence:
Cbaseline = m*CmaxIndex+b; % from regression
P = Cmax - Cbaseline; % normalize with baseline value

% Convert CmaxIndex to F0:
CmaxIndex = CmaxIndex+tRange(1)-1; % index into C
F0 = fs/(CmaxIndex-1); % Hz, convert index to lag with '-1'

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

% Bye!
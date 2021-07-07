function [th,f] = gammaChirpTheta(f,b,fr)
% This function calculates the value of theta used in the gammachirp function.
% Input [default]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate gammatone filter [0:1000]
%    b -- real scalar, bandwidth term of gammatone function [1]
%    fr -- real scalar, Hz, center frequency of gammatone function [100]
% Output:
%    th -- Nx1 real vector, rad, theta value in range [-pi/2,pi/2]
%    f -- Nx1 real vector, Hz, same as input

% Check inputs:
if nargin<3
   fr = 100;
end;
if nargin<2
   b = 1;
end;
if nargin<1
   f = [0:1000]';
end;

% Calculate ERB frequency of fr:
ERB_fr = f2ERB(fr);
ERB_fr = max(f2ERB(0),ERB_fr); % ensure positive ERB_fr (can be < 0 if Pgcp low enough)

% Calculate numerator and denominator:
num = f-fr;
den = b*ERB_fr;

% Calculate theta:
th = atan(num/den);

return;

% Bye!

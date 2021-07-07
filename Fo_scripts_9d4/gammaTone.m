function [GT,f,parametersGamma] = gammaTone(f,parametersGamma,checkParametersFlag)
% This function creates a gammatone filter magnitude response.
% Input [default values]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate gammatone filter [0:1000]
%         Note: if f is empty, the function adds defaults to parametersGamma and returns.
%    parametersGamma -- struct of gammachirp parameters (defaults from gammaChirpCompress.m):
%       .n1 -- real scalar, order of gammatone function
%       .b1 -- real scalar, bandwidth term of gammatone function
%       .fr1 -- real scalar, Hz, center frequency of gammatone function
%    checkParametersFlag -- logical scalar, [true]=check parameters, false=don't (faster)
% Output:
%    GT -- Nx1 real vector, magnitude spectrum of gammatone filter at frequencies f.
%    f -- Nx1 real vector, Hz, same as input.
%
% Note: GT at center frequency = unity

% Mark Skowronski, December 3, 2012

% Get defaults, init outputs (ignore f which is also input):
[GCC,GC,junk,fp1,GT,GCpeak,defaults] = gammaChirpCompress([]);

% Check inputs:
if nargin<3
   checkParametersFlag = true;
end;
if nargin<2
   parametersGamma = defaults;
elseif checkParametersFlag
   % Set default parameters to missing fields:
   parametersGamma = checkParameters(parametersGamma,defaults);
end;
if nargin<1
   f = [0:1000]'; % Hz
elseif isempty(f),
   % No frequencies to analyze, so return default parameters:
   return;
end;

% Ensure f is COLUMN vector:
f = f(:);

% Get ERB(fr1):
ERB_fr = f2ERB(parametersGamma.fr1); % ERB units (=Hz)

% Calculate bPrime:
bPrime = 2*pi*parametersGamma.b1*ERB_fr;

% Calculate numerator and denominator:
num = bPrime^parametersGamma.n1;
den = (bPrime^2 + (2*pi*(f-parametersGamma.fr1)).^2).^(parametersGamma.n1/2);

% Calcualte GT:
GT = num./den;

return;

% Bye!
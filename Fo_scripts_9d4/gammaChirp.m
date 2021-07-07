function [GC,f,fp1,GT,GCpeak,parametersGamma] = gammaChirp(f,parametersGamma,checkParametersFlag)
% This function creates a gammachirp filter magnitude response.
% Input [default values]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate gammachirp filter [0:1000]
%         Note: if f is empty, the function adds defaults to parametersGamma and returns.
%    parametersGamma -- struct of gammachirp parameters (defaults from gammaChirpCompress.m):
%       .n1 -- real scalar, order of gammatone function
%       .b1 -- real scalar, bandwidth term of gammatone function
%       .c1 -- real scalar, chirp factor
%       .fr1 -- real scalar, Hz, center frequency of gammatone function
%    checkParametersFlag -- logical scalar, [true]=check parameters, false=don't (faster)
% Output:
%    GC -- Nx1 real vector, magnitude spectrum of gammachirp filter at frequencies f.
%    f -- Nx1 real vector, Hz, same as input.
%    fp1 -- real scalar, Hz, frequency of GC peak
%    GT -- Nx1 real vector, magnitude spectrum of gammatone filter at frequencies f.
%    GCpeak -- real scalar, peak of GC before normalization (GT peak = unity)
%    parametersGamma -- same as input, with missing field names added from defaults
%
% Note: GC at peak frequency = unity, and peak frequency different from fr1 when c1 not equal to zero.
% The gammachirp function equals the gammatone function when c1 = 0.

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

% Get gammatone response:
GT = gammaTone(f,parametersGamma,false);

% Check for zero chirp:
if parametersGamma.c1==0,
   GC = GT;
   fp1 = parametersGamma.fr1;
   GCpeak = 1;
   return;
end;

% Calculate theta term, stage 1:
th = gammaChirpTheta(f,parametersGamma.b1,parametersGamma.fr1);

% Calculate peak value of GC (GT peak at unity):
num = exp(parametersGamma.c1*atan(parametersGamma.c1/parametersGamma.n1));
den = (1+(parametersGamma.c1/parametersGamma.n1)^2)^(parametersGamma.n1/2);
GCpeak = num/den;

% Calculate normalized gammachirp function:
GC = GT.*exp(parametersGamma.c1*th)/GCpeak;

% Calculate gammachirp peak frequency:
ERB_fr = f2ERB(parametersGamma.fr1);
fp1 = parametersGamma.fr1+parametersGamma.c1*parametersGamma.b1*ERB_fr/parametersGamma.n1;

return;

% Bye!
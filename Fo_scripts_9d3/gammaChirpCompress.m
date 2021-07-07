function [GCC,GC,f,fp1,GT,GCpeak,parametersGamma] = gammaChirpCompress(f,parametersGamma,...
   checkParametersFlag)
% This function creates a compressed gammachirp filter magnitude response, which is a cascade of a
% gammachirp function (stage 1) and a high-pass asymmetric function (stage 2).
% Input [default values]:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate gammatone filter [0:1000].
%         Note: if f is empty, the function adds defaults to parametersGamma and returns.
%    parametersGamma -- struct of parameters for compressed gammachirp function
%       Note: if parametersGamma is empty or not input, parametersGamma returned with default values.
%       Any missing fields are set to default values.
%       .auditoryFilterType -- string, ['GCC'],'roex' of auditory filter method to use
%       .n1 -- real scalar, order of gammatone function [4]
%       .b1,.b2 -- real scalar, bandwidth term of gammatone function [1.81,2.17]
%       .fr1 -- real scalar, Hz, center frequency of gammatone function [100]
%       .c1,.c2 -- real scalar, chirp factor [-2.96,2.20]
%       Note: if c2 is zero, GCC = GC and the function returns
%       .frat0,.frat1 -- real scalar, frequency ratio constant and slope, relating stage 2 asymptotic
%                   frequency to stage 1 output level, Pgcp [0.466,0.0109]
%       .Pgcp -- real scalar, dB, output power from passive gammachirp filter [50]
%    checkParametersFlag -- logical scalar, [true]=check parameters, false=don't (faster)
% Output:
%    GCC -- Nx1 real vector, magnitude spectrum of compressed gammachirp filter at frequencies f.
%    GC -- Nx1 real vector, magnitude spectrum of gammachirp filter at frequencies f.
%    f -- Nx1 real vector, Hz, same as input.
%    fp1 -- real scalar, Hz, frequency of GC peak
%    GT -- Nx1 real vector, magnitude spectrum of gammatone filter at frequencies f.
%    GCpeak -- real scalar, peak of GC before normalization (GT peak = unity)
%    parametersGamma -- same as input, with missing field names added from defaults
%
% Note: GC at peak frequency = unity, and peak frequency different from fr when c not equal to zero.
% The gammachirp function equals the gammatone function when c = 0.

% Mark Skowronski, December 3, 2012

% Set defaults:
defaults(1).auditoryFilterType = 'GCC';
defaults(1).n1 = 4; % gamma function order, stage 1
defaults(1).b1 = 1.81; % gammachirp bandwidth factor, stage 1
defaults(1).c1 = -2.96; % gammachirp factor, stage 1
defaults(1).fr1 = 100; % Hz, gammatone center frequency, stage 1
defaults(1).b2 = 2.17; % gammachirp bandwidth factor, stage 2
defaults(1).c2 = 2.20; % gammachirp factor, stage 2
defaults(1).frat0 = 0.466; % HP-AF frequency ratio constant, stage 2
defaults(1).frat1 = 0.0109; % HP-AF frequency ratio slope, stage 2
defaults(1).Pgcp = 50; % dB, Power output from normalized passive gammachirp filter

% Init outputs (ignore f and parametersGamma that are inputs):
GCC = [];
GC = [];
fp1 = [];
GT = [];
GCpeak = [];

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

% Get gammatone and passive gammachirp responses and peak frequency of stage 1:
[GC,f,fp1,GT,GCpeak] = gammaChirp(f,parametersGamma,false);

% Check for zero chirp in stage 2:
if parametersGamma.c2==0,
   GCC = GC;
   return;
end;

% Calculate stage 2 magnitude spectrum:
HPAF = calculateHPAF(f,fp1,parametersGamma);

% Combine GC and HPAF:
GCC = GC.*HPAF;

return;

% Bye!
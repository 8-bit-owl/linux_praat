function [excitationOut,FB,parameters] = signal2excitation(x,fs,fcAll,parameters,FB,checkFlag)
% This function calculates an excitation pattern for input x at frequencies f using a bank of
% auditory filters using one of the following methods:
% --Rounded exponential (roex), Glasberg & Moore 1990
% --Compressive gammachirp (GCC), Patterson, Unoki, & Irino 2003
% Each filter is divided into 2 stages:
% --Stage 1:
%    roex: symmetric roex function (input power assumed 51 dB/ERBrate)
%    GCC: passive gammachirp function
% --Stage 2:
%    roex: roex function, parameters calculated from power calculated from stage 1
%    GCC: highpass-asymmetric function, parameters calculated from power calculated from stage 1
% Note: auditory filter models are constrained such that filter bandwidth is not too wide nor too
% narrow for extreme input levels.
% Outer and middle ear transfer functions are applied to the signal spectrum before excitation
% calculation.
%
% Input [defaults]:
%    x -- Lx1 real vector, time domain signal
%    fs -- real scalar, Hz, sampling rate of x
%    fcAll -- Mx1 real vector, Hz, center frequencies of auditory filters [100]
%    parameters -- struct of parameters for auditory filter models
%       .auditoryFilterType -- string, ['GCC'],'roex' of auditory filter method to use
%       .FFTsize -- integer scalar, FFT size of x (if FFTsize<length(x), FFTsize set to [length(x)])
%       .windowType -- string of FFT window name: ['hamming'],'hanning','blackman','rectangular'
%       .transducer -- string, name of acoustic transducer, ['freefield'],'er2', see
%                      outerEarConversion.m
%       Note: Defaults for compressed gammachirp filter from gammaChirpCompress.m.
%       .n1 -- real scalar, order of gammatone function
%       .b1,.b2 -- real scalar, bandwidth term of gammatone function
%       .fr1 -- real scalar, Hz, center frequency of gammatone function (set by values in fcAll)
%       .c1,.c2 -- real scalar, chirp factor
%       .frat0,.frat1 -- real scalar, frequency ratio constant and slope, relating stage 2
%                        asymptotic frequency to stage 1 output level, Pgcp
%       .Pgcp -- real scalar, dB, output power from passive gammachirp filter
%    FB -- NxM real matrix, auditory filter bank, N bin frequencies, M filterbank center
%          frequencies.  Each filter is a magnitude squared spectrum. Set FB=[] to design FB from
%          parameters.
%    checkFlag -- logical scalar, [true]=check parameters, false=don't (faster)
% Output:
%    excitationOut -- Mx1 real vector, power units (not dB), auditory filterbank output
%    FB -- same as input, or designed from parameters and x spectrum if empty. Note: roex filter is
%          specified as squared magnitude function in Glasberg and Moore reference, while GCC filter
%          is specified as magnitude function (not squared) in Patterson et al. reference.
%    parameters -- same as input with parameters compared to defaults if flagged, or set to
%                  defaults if empty or not input.
%
% References: 
% Glasberg, B. R. and B. C. J. Moore, "Derivation of auditory filter shapes from notched-noise
% data," Hearing Research, vol. 47, pp. 103-138, 1990
% Patterson, R. D., M. Unoki, and T. Irino, "Extending the domain of center frequencies for the
% compressive gammachirp auditory filter," JASA, vol. 114(3), pp. 1529-1542, 2003

% Mark Skowronski, February 13, 2013

% Parameters:
minExcitation = 1e-10; % power units

% Get GCC defaults, add roex and FFT parameters:
defaults = getDefaults();
defaults(1).auditoryFilterType = 'GCC';
defaults(1).FFTsize = length(x);
defaults(1).windowType = 'hamming';
defaults(1).transducer = 'freefield';

% Check inputs:
if nargin<6
   checkFlag = true;
end;
if nargin<5
   FB = []; % empty, so design FB from parameters and x spectrum
end;
if nargin<4
   parameters = defaults;
elseif checkFlag || isempty(parameters)
   % Set default parameters to missing fields:
   parameters = checkParameters(parameters,defaults);

   % Check FFTsize:
   if parameters.FFTsize < length(x),
      parameters.FFTsize = length(x);
   end;
end;
if nargin<3
   fcAll = 100; % Hz
end;
if nargin<2
   % Cannot continue without x and fs input, so report error and return:
   error('Inputs x and fs required.');
end;

% Ensure inputs are COLUMNS:
x = x(:);
fcAll = fcAll(:);

% Construct power spectrum of x (COLUMN vector):
[X,fX,N] = getSpectrum(x,fs,parameters.FFTsize,parameters.windowType);

% Apply outer and middle ear transfer functions to X:
X = outerEarConversion(X,fX,parameters.transducer);
X = middleEarConversion(X,fX);

% Design FB if necessary:
if isempty(FB),
   FB = designFB(parameters,fcAll,fX,X,N);
end;

% Calculate excitation pattern using GCCbank:
excitationOut = integratePower(X,FB,N); % power units

% Add minimum excitation power:
excitationOut = excitationOut + minExcitation;

return;



function [GCbank,fp1All] = makeGCbank(fr1All,fX,parametersGamma)
% This function creates a matrix of passive gammachirp filter magnitudes.
% Input:
%    fr1All -- Mx1 real vector, Hz, gammatone center frequencies
%    fX -- Nx1 real vector, Hz, N bin frequencies of magnitude spectrum
%    parametersGamma -- struct of parameters for constructing gammachirp filters
%    Note: .fr1 is replaced with values in fr1All for each filter in the bank
% Output:
%    GCbank -- NxM real matrix, gammatone filter magnitude, N bin frequencies, M filterbank center
%              frequencies
%    fp1All -- Mx1 real vector, Hz, peak frequencies of passive gammachirp filters

% Init output:
M = length(fr1All);
N = length(fX);
GCbank = zeros(N,M);
fp1All = zeros(M,1);

% Construct each passive gammachirp filter, save fp1:
for p=1:M, % for each gammatone center frequency
   % Set gammatone center frequency in parametersGamma:
   parametersGamma.fr1 = fr1All(p);

   % Get passive gammachirp filter:
   [GCbank(:,p),fX,fp1All(p)] = gammaChirp(fX,parametersGamma,false);
end;

return;


function PgcpAll = zeroHLtruncate(PgcpAll,fp1All)
% This function truncates the stage 1 output in PgcpAll to the output of a 0 dB HL tone at
% frequencies fp1All.

% Get dB SPL level for 0 dB HL at fp1All:
dBSPL = zeroHL2SPL(fp1All,'ER3A'); % use ER3A transducer

% Truncate PgcpAll at dBSPL:
PgcpAll = max(PgcpAll,dBSPL);

return;


function [GCCbank,HPAFbank] = makeGCCbank(GCbank,fp1All,fX,PgcpAll,parametersGamma)
% This function constructs the stage 2 HP-AF filters and combines them with GCbank.
% Input:
%    GCbank -- NxM real matrix, gammatone filter magnitude, N bin frequencies, M filterbank center
%              frequencies
%    fp1All -- 1xM real vector, Hz, passive gammachirp filter peak frequency
%    fX -- Nx1 real vector, Hz, N bin frequencies of magnitude spectrum
%    PgcpAll -- 1xM real vector, dB, output power from each filter in the bank
%    parametersGamma -- struct of compressive gammachirp filter parameters (see gammChirpCompress.m)
% Output:
%    GCCbank -- same size as GCbank, includes stage 2 magnitude response

% Init output:
GCCbank = zeros(size(GCbank));
HPAFbank = zeros(size(GCbank));

% Construct stage 2 for each filter, combine with GC:
for p=1:length(fp1All),
   % Store Pgcp in parametersGamma:
   parametersGamma.Pgcp = PgcpAll(p);

   % Calculate HP-AF magnitude spectrum:
   HPAF = calculateHPAF(fX,fp1All(p),parametersGamma);
   HPAFbank(:,p) = HPAF;

   % Combine HPAF and GC:
   GCCbank(:,p) = GCbank(:,p).*HPAF;
end;

return;


function FB = makeRoexBank(fcAll,fX,PEstimate)
% This function calculates a bank of roex filters.
% Input:
%    fcAll -- Mx1 real vector, Hz, center frequencies of auditory filters
%    fX -- Nx1 real vector, Hz, N bin frequencies of magnitude spectrum
%    PEstimate -- Mx1 real vector or real scalar, dB, equivalent input noise level
% Output:
%    FB -- NxM real matrix, squared magnitude spectra, N bin frequencies, M filterbank center
%          frequencies

% Parameters:
pMin = 3; % minimum p value (sets maximum filter bandwidth)
pMax = 40; % maximum p value (sets minimum filter bandwidth)

% Ensure inputs are COLUMN vectors:
fcAll = fcAll(:);
fX = fX(:);
PEstimate = PEstimate(:);

% Calculate ERB for each fc:
ERBAll = f2ERB(fcAll); % Hz

% Calculate p for upper tail, constrain:
puAll = 4*fcAll./ERBAll;
puAll = max(pMin,min(pMax,puAll)); % constrain puAll to be in range [pMin,pMax]

% Calculate p for lower tail, constrain:
pu1k = 4*1000/f2ERB(1000); % pu at 1000 Hz
plAll = puAll-0.38*puAll/pu1k.*(PEstimate-51);
plAll = max(pMin,min(pMax,plAll)); % constrain plAll to be in range [pMin,pMax]

% Init output:
FB = zeros(length(fX),length(fcAll));

% Design filter for each fc:
for p=1:length(fcAll),
   FB(:,p) = roexFunction(fX,fcAll(p),puAll(p),plAll(p));
end;

return;


function FB = designFB(parameters,fcAll,fX,X,N)
% This function designs the stage 1 and stage 2 filterbanks.  The output is the squared magnitude
% spectrum of each filter, one filter per COLUMN.

% Construct stage 1 filter bank:
switch parameters.auditoryFilterType
   case 'roex', % standard roex filter, 51 dB/ERB power
      FB = makeRoexBank(fcAll,fX,51); % squared magnitude spectrum
      fp1All = fcAll; % stage 1 peak frequency = fc for roex
   otherwise, % passive gammachirp filter
      [GCbank,fp1All] = makeGCbank(fcAll,fX,parameters);
      FB = GCbank.^2; % convert to squared magnitude
end;

% Calculate stage 1 output power:
PEstimate = integratePower(X,FB,N); % power units
PEstimate = 10*log10(PEstimate); % dB

% Limit stage 1 output power to 0 dB HL at fp1All:
dBSPL = zeroHL2SPL(fp1All);
PEstimate = max(dBSPL,PEstimate);

% Construct stage 2 filter bank:
switch parameters.auditoryFilterType
   case 'roex', % roex filter, parameters modified by PEstimate
      FB = makeRoexBank(fcAll,fX,PEstimate); % squared magnitude spectrum
   otherwise, % stage 1 filterbank + stage 2 HP-AF
      GCCbank = makeGCCbank(GCbank,fp1All,fX,PEstimate,parameters);
      FB = GCCbank.^2; % convert to squared magnitude
end;

return;

% Bye!

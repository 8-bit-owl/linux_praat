function [Ls,Es,Ln,En,fcAll,parameters,FBa,fX] = loudnessChen(xSignal,xNoise,fs,parameters)
% This function calculates excitation functions and total loudness for a signal in masking noise
% according to the method in [1]. The signal and noise input are converted to power spectra and
% excitation functions using a parallel roex filter architecture: 1) passive wideband roex, output
% power sets gain for active filter, 2) active narrowband roex, gain decreases as passive output
% power increases (compression).  Because the excitation function includes compression, total
% loudness is calculated directly from the excitation function instead of indirectly through the
% calculation of a specific loudness function as in previous studies from the Moore Lab.
%
% Notes:
% ---Set xNoise to 0 or [] to calculate noise-free loudness and excitation.
%
% Input:
%    xSignal,xNoise -- Lx1 real vectors, time domain data for signal and noise
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of design parameters [defaults]
%       .fcRange -- 1x2 real vector, Hz, [min max] center frequencies [40 17000]
%       .ERBspacing -- real scalar, ERBrate, spacing of center frequencies in ERBrate frequency.
%                      Note: the first filter center frequency is fcRange(1), and the last filter
%                      center frequency is the largest multiple of ERBspacing above fcRange(1) that
%                      is less than or equal to fcRange(2). [0.1]
%       .FFTsize -- integer scalar, FFT size of x [max(2^14,length(x))]
%       .windowType -- string of FFT window name: ['hamming'],'hanning','blackman','rectangular'
%       .transducer -- string, name of acoustic transducer, ['freefield'],'er2', see
%                      outerEarConversion.m
% Output:
%    Ls,Ln -- real scalars, sones, loudness for signal, noise
%    Es,En -- Mx1 real vectors, power units, excitation functions for signal, noise
%    fcAll -- Mx1 real vector, Hz, center frequencies of M auditory filters, design from parameters
%    parameters -- same as input with parameters compared to defaults if flagged, or set to
%                  defaults if empty or not input.
%
% References:
% [1] Chen, Z., G. Hu, B. R. Glasberg, and B. C. J. Moore, "A new method of calculating auditory
% excitation patterns and loudness for steady sounds," Hearing Research, vol. 282, pp. 204-215, 2011

% Mark Skowronski, February 21, 2013

% Parameters:
minLoudness = 1e-5; % sones

% Set defaults:
defaults(1).fcRange = [40 17000]; % Hz
defaults(1).ERBspacing = 0.1; % ERBrate
defaults(1).FFTsize = max(2^14,length(xSignal));
defaults(1).windowType = 'hamming';
defaults(1).transducer = 'freefield';

% Check inputs:
if nargin<4
   parameters = defaults;
end;
if nargin<3
   error('Must input xSignal, xNoise, and fs.');
elseif sum(abs(xNoise))==0, % xNoise is 0 or []
   xNoise = zeros(size(xSignal)); % convert xNoise to vector of zeros
elseif length(xSignal(:))~=length(xNoise(:)),
   error('xSignal and xNoise must be same length.');
end;

% Check that FFTsize is >= length(xSignal),
if parameters.FFTsize < length(xSignal),
   parameters.FFTsize = length(xSignal);
end;

% Make fcAll:
fcAll = makeFc(parameters.fcRange,parameters.ERBspacing);

% Construct signal, noise, and total spectra:
[Xsignal,fX,N] = getSpectrum(xSignal,fs,parameters.FFTsize,parameters.windowType);
Xnoise = getSpectrum(xNoise,fs,parameters.FFTsize,parameters.windowType); % same fX and N from signal
%Xtotal = getSpectrum(xSignal+xNoise,fs,parameters.FFTsize,parameters.windowType); % used for active FB design

% Apply outer ear conversion:
Xsignal = outerEarConversion(Xsignal,fX,parameters.transducer);
Xnoise = outerEarConversion(Xnoise,fX,parameters.transducer);
%Xtotal = outerEarConversion(Xtotal,fX,parameters.transducer);

% Apply middle ear conversion:
Xsignal = middleEarConversion(Xsignal,fX);
Xnoise = middleEarConversion(Xnoise,fX);
%Xtotal = middleEarConversion(Xtotal,fX);

% Construct passive filter bank:
FBp = makeFBp(fcAll,fX);

% Calculate passive filter output excitation:
EpSignal = integratePower(Xsignal,FBp,N);
EpNoise = integratePower(Xnoise,FBp,N);
%EpTotal = integratePower(Xtotal,FBp,N); % used for active FB design
EpTotal = EpSignal + EpNoise;

% Construct active filter bank:
FBa = makeFBa(fcAll,fX,EpTotal);

% Calculate active filter output excitation:
EaSignal = integratePower(Xsignal,FBa,N);
EaNoise = integratePower(Xnoise,FBa,N);
%EaTotal = integratePower(Xtotal,FBa,N);
EaTotal = EaSignal + EaNoise;

% Combine passive and active filter excitations:
Es = EpSignal+EaSignal;
En = EpNoise+EaNoise;
Et = EpTotal+EaTotal;

% Calculate loudness:
Lt = getLoudness(Et,0,fcAll);
Ls = getLoudness(Es,En,fcAll);
Ln = max(0,Lt-Ls); % ensure non-negative

return;


function FB = makeFBp(fcAll,fX)
% This function constructs the passive filter bank from Chen et al. [1].  Parameters for the
% equations relating filter bank parameters to fcAll are hard-coded.
% Input:
%    fcAll -- Mx1 real vector, Hz, center frequencies of auditory filters
%    fX -- Qx1 real vector, Hz, FFT bin frequencies
% Output:
%    FB -- QxM real matrix, power units, filter bank, one roex filter per COLUMN

% Get lengths of inputs:
M = length(fcAll); % number of center frequencies
Q = length(fX); % number of FFT bins in spectrum

% Init output:
FB = zeros(Q,M);

% Calculate filter bank parameters from fcAll:
tl = fcAll./(0.108*fcAll+2.33); % lower skirt bandwidth term, tail roex
tu = 15.6; % upper skirt bandwidth term, tail roex

% Construct filter for each fc:
for p=1:M,
   FB(:,p) = roexFunction(fX,fcAll(p),tu,tl(p));
end;

return;


function FB = makeFBa(fcAll,fX,E)
% This function constructs the active filter bank from Chen et al. [1].  Parameters for the
% equations relating filter bank parameters to fcAll are hard-coded.
% Input:
%    fcAll -- Mx1 real vector, Hz, auditory filter center frequencies, M filters
%    fX -- Qx1 real vector, Hz, FFT bin frequencies, Q bins
%    E -- Mx1 real vector, power units, excitation output from passive filters
% Output:
%    FB -- QxM real matrix, power units, filter bank, one roex filter per COLUMN

% Calculate filter bank frequencies:
Q = length(fX);
M = length(fcAll);

% Init output:
FB = zeros(Q,M);

% Calculate filter bank parameters from fcAll:
GdBmax = fcAll./(0.0191*fcAll+1.10);
pl = fcAll./(0.027*fcAll+5.44);
pu = 27.9;

% Convert E to dB:
EdB = 10*log10(E);

% Calculate active filter gain for EdB above and below 30 dB:
GdB = GdBmax.*(1 - 1./(1+exp(-.05*(EdB-(100-GdBmax)))) + 1./(1+exp(.05*(100-GdBmax)))); % below 30 dB
t = find(EdB>30);
if ~isempty(t),
   GdB(t) = GdB(t)-0.003*(EdB(t)-30).^2; % above 30 dB
end;

% Convert EdB to power units:
G = 10.^(GdB/10);

% Construct filter for each fc:
for p=1:M,
   W = roexFunction(fX,fcAll(p),pu,pl(p)); % unity gain
   FB(:,p) = G(p)*W; % scale by G
end;

return;


function L = getLoudness(Es,En,fcAll)
% This function converts excitation to loudness according to Chen et al. [1].  Includes hard-coded
% parameters.

% Calculate filter bank parameters from fcAll:
C = 1.53e-8; % conversion factor to sones
Ef = f2ERBrate(fcAll)/f2ERBrate(1000)-1; % normalized ERBrate frequency
KdB = 6.51*Ef.^2-1.93; % dB, SNR at threshold in noise
K = 10.^(KdB/10); % power units

% Combine excitations:
E = max(0,Es-K.*En); % truncate at zero

% Integrate:
df = diff(f2ERBrate(fcAll)); % ERBrate
df(end+1) = df(end); % repeat last value to match vector lengths
L = C*sum(E.*df);

return;


% Bye!
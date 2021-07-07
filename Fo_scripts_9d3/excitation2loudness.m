function [Lsig,Lnoise,sLsig,sLnoise] = excitation2loudness(Esig,Enoise,fcAll)
% This function converts signal and noise excitations, Esig and Enoise, into specific loudness
% functions and total loudness.  Set Enoise=0 for noise-free loudness calculation of Esig.
% Input:
%    Esig -- Mx1 real vector, power units, excitation function for signal of interest, M auditory
%            filters
%    Enoise -- Mx1 real vector or scalar, power units, excitation function for masking noise
%    fcAll -- Mx1 real vector, Hz, center frequencies of auditory filters
% Output:
%    Lsig -- real scalar, sone, loudness of signal of interest, integrated from specific loudness
%    Lnoise -- real scalar, sone, loudness of noise, integrated from specific loudness
%    sLsig -- Mx1 real vector, sone/ERBrate, specific loudness of signal of interest
%    sLnoise -- Mx1 real vector, sone/ERBrate, specific loudness of noise
%
% Reference: Moore, Glasberg, and Baer, "A model for prediction of thresholds, loudness, and partial
% loudness," JAES 1997

% Mark Skowronski, February 12, 2013

% Ensure inputs are COLUMN vectors:
Esig = Esig(:);
Enoise = Enoise(:);
fcAll = fcAll(:);
M = length(Esig); % number of auditory filters used to construct excitation function
I = ones(M,1); % unity vector
Enoise = Enoise.*I; % vectorize scalar Enoise

% Parameters (constants above 500 Hz):
C = 0.047; % scale factor such that loudness of 1-kHz tone presented binaurally in free field at 40 dB SPL = 1 sone
Ethrq = 2.31; % power units, excitation level of tone at threshold in quiet (no noise)
G = 1; % unitless, relative gain of excitation due to cochlear amplifier
alphaExp = 0.2; % unitless, excitation exponent
A = 2*Ethrq; % power units, excitation additive factor
KdB = -3; % dB, SNR scale factor, Ethrn = K*Enoise+Ethrq, where K is in power units
K = 10^(KdB/10); % convert from dB scale to power units scale

% Set data from Figures for low-frequency variables:
EthrqYPlot = [26.2,20.2,14.5,6.3,10*log10(Ethrq)]; % dB, excitation at threshold in quiet (no noise), Figs. 4 and 8
EthrqXPlot = [52,74,108,253,500]; % Hz, Figs. 4 and 8, note: log frequency scale
alphaExpYPlot = [.2669,.2504,.2360,.2216,.2101,.2]; % unitless, Fig. 6
alphaExpXPlot = [-25,-20,-15,-10,-5,0]; % dB, G value, Fig. 6
AYPlot = [8.775,7.529,6.574,5.744,5.121,2*Ethrq]; % power units, Fig. 7
AXPlot = [-25,-20,-15,-10,-5,0]; % dB, G value, Fig. 7
KdBYPlot = [13.41,7.14,4.79,0.310,-1.069,-2.448,KdB]; % dB, SNR term, Fig. 9
KdBXPlot = [50,80,100,200,300,500,900]; % Hz, frequency, Fig. 9

% Interpolate for Ethrq for low frequencies from fcAll, note: log frequency scale, so use log(fcAll).
lowFreqIndex = find(fcAll<500); % index into fcAll, E
EthrqLowFreqdB = interp1(log10(EthrqXPlot),EthrqYPlot,log10(fcAll(lowFreqIndex)),'spline','extrap');
EthrqLowFreq = 10.^(EthrqLowFreqdB/10); % power units

% Interpolate for KdB for low frequencies from fcAll, note: log frequency scale, so use log(fcAll).
lowFreqIndex900 = find(fcAll<900); % index into fcAll, E
KdBLowFreq = interp1(log10(KdBXPlot),KdBYPlot,log10(fcAll(lowFreqIndex900)),'spline','extrap');
KLowFreq = 10.^(KdBLowFreq/10); % power units

% Calculate G from Ethrq:
GLowFreq = Ethrq./EthrqLowFreq; % unitless
GLowFreqdB = 10*log10(GLowFreq); % dB

% Interpolate for alpha and A for low frequencies from G:
alphaExpLowFreq = interp1(alphaExpXPlot,alphaExpYPlot,GLowFreqdB,'spline','extrap');
ALowFreq = interp1(AXPlot,AYPlot,GLowFreqdB,'spline','extrap');

% Init signal specific loudness vectors:
sLsig = zeros(M,1); % sone/ERBrate

% Init model parameters as COLUMN vectors of length M:
EthrqAll = Ethrq*I;
GAll = G*I;
alphaExpAll = alphaExp*I;
AAll = A*I;
KAll = K*I;

% Include low-frequency values for model parameters:
EthrqAll(lowFreqIndex) = EthrqLowFreq;
GAll(lowFreqIndex) = GLowFreq;
alphaExpAll(lowFreqIndex) = alphaExpLowFreq;
AAll(lowFreqIndex) = ALowFreq;
KAll(lowFreqIndex900) = KLowFreq; % constant above 900 Hz, not 500 Hz like other parameters

% Calculate excitation level of tone at threshold in noise:
EthrnAll = KAll.*Enoise+EthrqAll; % power units, p. 230, first paragraph of reference

% Calculate total excitation:
Etotal = Esig+Enoise; % power units

% Calculate signal specific loudness: Esig+Enoise>=1e10, above threshold:
iha = find(Etotal>=1e10 & Esig>=EthrnAll); % index into Esig, Enoise, specific loudness functions
if ~isempty(iha),
   sLsig(iha) = equation19(Esig(iha),Enoise(iha),C,EthrnAll(iha),EthrqAll(iha),...
      KAll(iha),GAll(iha),AAll(iha),alphaExpAll(iha));
end;

% Calculate signal specific loudness: Esig+Enoise>=1e10, below threshold:
ihb = find(Etotal>=1e10 & Esig<EthrnAll); % index into Esig, Enoise, specific loudness functions
if ~isempty(ihb),
   sLsig(ihb) = equation20(Esig(ihb),Enoise(ihb),C,EthrnAll(ihb),EthrqAll(ihb),...
      KAll(ihb),GAll(ihb),AAll(ihb),alphaExpAll(ihb));
end;

% Calculate signal specific loudness: Esig+Enoise<1e10, above threshold:
ila = find(Etotal<1e10 & Esig>=EthrnAll);
if ~isempty(ila),
   sLsig(ila) = equation17(Esig(ila),Enoise(ila),C,EthrnAll(ila),EthrqAll(ila),...
      KAll(ila),GAll(ila),AAll(ila),alphaExpAll(ila));
end;

% Calculate signal specific loudness: Esig+Enoise<1e10, below threshold:
ilb = find(Etotal<1e10 & Esig<EthrnAll); % index into Esig, Enoise, specific loudness functions
if ~isempty(ilb),
   sLsig(ilb) = equation18(Esig(ilb),Enoise(ilb),C,EthrnAll(ilb),EthrqAll(ilb),...
      KAll(ilb),GAll(ilb),AAll(ilb),alphaExpAll(ilb));
end;

% Calculate noise specific loudness:
if sum(Enoise)==0, % noise-free case
   sLnoise = zeros(M,1); % sone/ERBrate
else
   % Calculate specific loudness for Etotal (recursive):
   [~,~,sLtotal] = excitation2loudness(Esig+Enoise,0,fcAll);
   sLnoise = sLtotal-sLsig;
end;

% Calculate total loudness as integral of specific loudness over ERBrate:
Lsig = specificLoudness2TotalLoudness(sLsig,fcAll);
Lnoise = specificLoudness2TotalLoudness(sLnoise,fcAll);

return;


function sLsig = equation17(Esig,Enoise,C,Ethrn,Ethrq,K,G,A,alphaExp)
% This function implements Eq. 17 which relates signal specific loudness to excitation for
% Etotal<1e10, above threshold.

part1 = C*((G.*(Esig+Enoise)+A).^alphaExp-A.^alphaExp);
part2 = (Ethrn./Esig).^0.3; % attenuation factor
part3 = C*((G.*(Enoise.*(K+1)+Ethrq)+A).^alphaExp-(G.*Ethrq+A).^alphaExp);

sLsig = part1-part2.*part3;

return;


function sLsig = equation18(Esig,Enoise,C,Ethrn,Ethrq,K,G,A,alphaExp)
% This function implements Eq. 18 which relates signal specific loudness to excitation for
% Etotal<1e10, below threshold.

part1 = C*((G.*(Esig+Enoise)+A).^alphaExp-(G.*Enoise+A).^alphaExp);
part2 = (2*Esig./(Esig+Ethrn)).^1.5;
part3 = ((G.*Ethrq+A).^alphaExp-A.^alphaExp)./((G.*(Enoise.*(K+1)+Ethrq)+A).^alphaExp-(G.*Enoise+A).^alphaExp);

sLsig = part1.*part2.*part3;

return;


function sLsig = equation19(Esig,Enoise,C,Ethrn,Ethrq,K,G,A,alphaExp)
% This function implements Eq. 19 which relates signal specific loudness to excitation for
% Etotal>=1e10, above threshold.

C2 = C/sqrt(1.04e6);

part1 = C2*sqrt(Esig+Enoise); % approximate loudness of total signal
part2 = (Ethrn./Esig).^0.3; % attenuation factor
part3 = C2*(sqrt(Enoise.*(K+1)+Ethrq)-(G.*Ethrq+A).^alphaExp+A.^alphaExp);

sLsig = part1-part2.*part3;

return;


function sLsig = equation20(Esig,Enoise,C,Ethrn,Ethrq,K,G,A,alphaExp)
% This function implements Eq. 20 which relates signal specific loudness to excitation for
% Etotal>=1e10, below threshold.

part1 = C*(sqrt(Esig+Enoise)-sqrt(Enoise));
part2 = (2*Esig./(Esig+Ethrn)).^1.5;
part3 = ((G.*Ethrq+A).^alphaExp-A.^alphaExp)./(sqrt(Enoise.*(K+1)+Ethrq)-sqrt(Enoise));

sLsig = part1.*part2.*part3;

return;

% Bye!
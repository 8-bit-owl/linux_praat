function HPAF = calculateHPAF(f,fp1,parametersGamma)
% This function calculates the high-pass asymmetric function (stage 2) of the compressive gammachirp 
% filter.
% Input:
%    f -- Nx1 real vector, Hz, frequencies at which to evaluate gammatone filter [0:1000].
%         Note: if f is empty, the function adds defaults to parametersGamma and returns.
%    fp1 -- real scalar, Hz, frequency of GC peak
%    parametersGamma -- struct of compressive gammachirp filter parameters (see gammChirpCompress.m)
% Output:
%    HPAF -- Nx1 real vector, magnitude spectrum of high-pass asymmetric function (level dependent)

% Compare Pgcp to 0 dB HL at fp1, clip if necessary:
dBSPL = zeroHL2SPL(fp1,'ER3A');
parametersGamma.Pgcp = max([parametersGamma.Pgcp,dBSPL]);

% Calculate level-dependent asymptotic frequency of stage 2:
fr2 = fp1*(parametersGamma.frat0+parametersGamma.frat1*parametersGamma.Pgcp); % Hz

% Calculate theta, stage 2:
th2 = gammaChirpTheta(f,parametersGamma.b2,fr2);

% Calculate HPAF:
HPAF = exp(parametersGamma.c2*th2);

return;

% Bye!
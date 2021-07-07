function [E_norm, ppd] = ppe2(F0)
% Calculates the (F0) pitch period entropy for a given input voice recording,
% using standard algorithm parameter values. Based on Praat's pitch
% extraction implementation.
%
% Usage:
% [E_norm, ppd] = ppe2(F0)
% Inputs
%    F0     - fundamental frequency of signal
%
% Outputs:
%    E_norm - pitch period entropy value
%    ppd    - pitch period density
%
% (c) 2008 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, E.J. Hunter, J. Spielman, L.O. Ramig (2008)
% Suitability of dysphonia measurements for telemonitoring of Parkinson’s
% disease
% IEEE Transactions on Biomedical Engineering
% 20180907- edits by Eric Hunter

% Obtain pitch series using Praat
% [t, F0] = praat_pitch(x, fs);

% Convert to perceptual pitch scale
pseries = 12*log2(F0/127.09);

% Whiten series and remove transient
if (length(pseries) > 20)
    a = arcov(pseries, 2);
    f = filter(a, 1, pseries);
    f = f(10:end);

    % Calculate spread measure
    xbins  = linspace(-4.3, 2.7, 60);
    ppd    = hist(f, xbins);
    ppd    = ppd/sum(ppd);
    E_norm = entropy(ppd)/log(length(ppd));
else
    warning('F0 series not long enough to calculate PPE');
    E_norm = NaN;
end

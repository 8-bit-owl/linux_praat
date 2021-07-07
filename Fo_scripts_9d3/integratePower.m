function P = integratePower(X,FB,N)
% This function combines the power spectrum in X and FB to calculate output power.
% Input:
%    X -- Qx1 real vector, power spectral density, lower half of spectrum
%    FB -- QxM real matrix, bank of filter squared magnitudes, Q bin frequencies, M filterbank
%          center frequencies
%    N -- integer scalar, FFT size
% Output:
%    P -- Mx1 real vector, power units (not dB), output power from each filter in the bank

% Calculate inner produce of X and FB based on odd or even N:
if mod(N,2)==0, % even length, X includes DC, lower spectrum, and Nyquist frequencies.
   P = calcIP(X,FB,1) + 2*calcIP(X,FB,[2:length(X)-1]) + calcIP(X,FB,length(X));
else % odd length, X includes DC and lower spectrum frequencies.
   P = calcIP(X,FB,1) + 2*calcIP(X,FB,[2:length(X)]);
end;
P = P/N; % scale by FFT size to convert P to average

return;


function P = calcIP(X,FB,t)
% This function calculates the inner product between X and FB for rows t.

P = FB(t,:)'*X(t);

return;


% Bye!
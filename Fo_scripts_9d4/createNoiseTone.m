function x = createNoiseTone(F0,fs,T,SNR,band)
% This function creates a tone at F0 of duration T seconds in the presence of bandlimited noise.
% Input:
%    F0 -- real scalar, Hz, frequency of tone
%    fs -- real scalar, Hz, sampling rate of x
%    T -- real scalar, sec, duration of x
%    SNR -- real scalar, dB, signal to noise ratio
%    band -- string, ['low'],'high','both', include noise below F0, above F0, or surrounding F0
% Output:
%    x -- Nx1 real vector, output signal, N = round(T*fs)

% Parameters:
filterOrder = 50; % Butterworth filter order

% Init output with tone at F0:
N = round(T*fs);
t = [0:N-1]'/fs; % sec, COLUMN vector
x = sin(2*pi*F0*t);

% Create noise filter:
switch band
   case 'both'
      % Band-stop filter:
      Wn = [0.90,1.1]*F0/(fs/2);
      [z,p,k] = butter(filterOrder,Wn,'stop');
   case 'high'
      % High-pass filter:
      Wn = 1.1*F0/(fs/2); % 10 percent above F0
      [z,p,k] = butter(filterOrder,Wn,'high');
   otherwise
      % Low-pass filter:
      Wn = .95*F0/(fs/2); % 10 percent below F0
      [z,p,k] = butter(filterOrder,Wn); % pole-zero-gain output
end;

% Create noise, filter, SNR scale factor:
xn = randn(size(x));
xn = xFilter(xn,z,p,k);
xPow = mean(x.^2);
xnPow = mean(xn.^2);
An = sqrt(xPow/xnPow/10^(SNR/10)); % noise amplitude scale factor

% Combine x and xn:
x = x+An*xn;

return;

% Bye!

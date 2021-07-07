function [ftable fm]  = octaves(N,b,base)
% OCTAVES calculates the IEC/ANSI standard octave band associated with the 
%   number N. 'b' is the 'bandwidth designator' or in other words the
%   fraction of octaves desired (use 1 for octaves, 3
%   for third octaves). The integer N determines the BAND DISTANCE FROM the
%   1000 HZ band (which is always N=30), which will vary for octaves and 
%   third octaves. For example, N=29 is one band BELOW the 1000 Hz band 
%   which would be the 500 Hz band for octaves, or the 794 Hz band for 
%   third octaves. When the bandwidth designator b=3, N is known as the
%   'band number'.
%   The center frequencies are calculated according to 'base', which can be
%   2 or 10. (10 is preferred for the standard, but not required. 2 gives
%   traditional audiometric octaves.)
% 
% [ftable fm]  = octaves(N,b,base)
%   'ftable' is a length(N) x 4 matrix. The first column contains the band
%   number (N), the second is the center frequency, the third is the
%   lower cutoff frequency, and the fourth is the upper cutoff frequency.
%   ALL FREQUENCY VALUES ARE GIVEN IN KHZ.
%   'fm' is the center (or midband) frequency in kHz.
% 
% The function is based on the IEC 1260/ANSI S1.11 standard.
% 
% Example:
%       octanalysis = octaves(-3:4,1);
%           (calculates octave bands from 125 Hz to 16 kHz)
%       oct3analysis = octaves(-10:13,3);
%           (calculates third-octave bands from 100 Hz to 20.2 kHz)

switch base
    case 10
        G = 10^(3/10); % for base 10 calculation
    case 2
        G = 2; % for base 2 calculation
end

% fr = 1000;

% according to IEC 1260:1995
% if mod(b,2)
%     fm = (G.^(N/b));
% else
%     fm = (G.^((2*N+1)/(2*b)));
% end

% according to ANSI S1.11-2004 (R2009)
% if b = 3 then N is considered the 'band number'
if mod(b,2)
    fm = (G.^((N-30)/b));
else
    fm = (G.^((2*N-59)/(2*b)));
end

f1 = (G^(-1/(2*b))).*fm;
f2 = (G^(1/(2*b))).*fm;

ftable(:,1) = N;
ftable(:,2) = fm;

ftable(:,3) = f1;
ftable(:,4) = f2;

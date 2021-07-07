function [siga,sig,t2,Fs,fname]=LoadNewWave(fname,Fs,time_steps)
% LoadNewWave gets a wave file, filters and returns
%  Inputs:
%     Fs-   sample rate of the arrays
%     fname- name of the wave file loaded
%     time_steps- array of two numbers identifying the window of time (in 
%     minutes)to read from the wave file. for example =[1 1.01] loads from 
%     minute 1 to minute 1.01.  
%  Returns:
%     siga- the raw signal from the wave file
%     sig-  filtered, normalized, zeroed signal for viewing, etc
%     t1-   time array corresponding to the above two signals (in seconds)
%     Fs-   sample rate of the arrays
%     fname- name of the wave file loaded
%
% Eric Hunter 20050615, adjusted 20060417.
if nargin<1 
   [fname]=uigetfile('*.wav;*.WAV','Open the sound file');
   [Y,Fs,nbits] = wavread(fname,100);
   time_steps=[0 0.3];
end
N1 = floor(time_steps(1)*Fs*60)+1;
N2 = floor(time_steps(2)*Fs*60);
[siga,Fs,nbits] = wavread(fname,[N1 N2-1]);
t2=(N1-1:N2-2)/(Fs);      % get time array in seconds

% arry = find(t1<20);
% t2 = t1(arry);
% sig_short = siga(arry);
% sig_short = sig_short(1:2:end);
% t2 = t2(1:2:end); 

try
%% filter the signal and send both filtered and raw back
fu = 2000;  %upper bandwidth to filter the signal
% sig = siga/max(abs(siga));  % normalize
% sig = sig-mean(sig);        % zero
c = fir1(250,[200]/(Fs/2),'high');  % filter
sig = filtfilt(c,1,siga);

catch
   sig=siga; 
end


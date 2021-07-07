function [sig_RMS,sig_dB,sig_fftE,t2,gain_volt,frame_sz]=GetRMSandDBandSpecEnergy(...
   siga,frm_sz_ms,frm_step_ms,Fs,step,gain_volt)
%GetRMSandDB takes the signal and finds the RMS and dB from it
%
%  requires:
%     siga-       array of the audio signal
%     frm_sz_ms-  frame size (in msec) of the window to get RMS over
%                 20ms -> 50Hz, 15-> 66.6Hz
%     frm_step_ms-frame step (in msec) 
%     Fs-         the sampling rate of the signal
%     step-       needed to keep the time index up to date
%   optional:
%     gain_volt-  voltage gain for calibrating to dB SPL, otherwise some
%                 default is used which means that dB is only relative
%
%  Returns:
%     sig_RMS-    the RMS of siga at steps of frm_sz_ms
%     sig_dB-     the dB of siga  at steps of frm_sz_ms
%     sig_fftE-   the total engergy in the FFT of the signal in same steps
%     t2-         the time array for both of the above in minutes
%     gain_volt-  the voltage gain used to get dB.  this is sent back in
%                 the instance that no gain was sent, then the user knows
%                 the value used.  a '0.5' is always at the end of this to
%                 indicate that an arbitrary value was used.
%
%  -Eric Hunter, 20060417

if nargin<6
   gain_volt=84692.5;
end

frame_sz=floor(Fs*(frm_sz_ms/1000)); 
frame_step = Fs*(frm_step_ms/1000);
sig_RMS=rms_eh(siga,frame_sz,frame_step);
sig_fftE=fftE_eh(siga,frame_sz,frame_step,Fs);

t1=(1:length(siga))/Fs;
% t2=t1(1:frame_sz:length(sig_RMS)*frame_sz);
t2 = step + frm_step_ms*(0:length(sig_RMS)-1)/(60*1000);
sig_RMS(sig_RMS <= 0) = .000000001;
sig_dB=20*log10(gain_volt.*sig_RMS);


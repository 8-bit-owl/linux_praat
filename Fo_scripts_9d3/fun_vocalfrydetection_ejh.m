function [time_vf,vf_segm]=fun_vocalfrydetection_ejh(filename,power_th,fignum)

% Copyright July 2016
% This program was developed by Carlos Ishi (ATR).
% This program detects vocal fry segments for an input wav file.
% The wav file name and a power threshold (in dB) have to be set. 
% The sampling rate of the wav file has to be 16 kHz.
% The power threshold has to be adapted for different recordings.
% Look at the power contour and power threshold (red line) after processing.
% The red line should be above noise power and below vocal fry pulse power.
% Other parameters related to vocal fry detection can also be modified below.
% The variables starttime and endtime retain the start and end time of the detected vocal fry segments.
% Reference: 
% Ishi, C.T., Sakakibara, K-I., Ishiguro, H., Hagita, N. (2008). A method for automatic detection of vocal fry. IEEE Transactions on Audio, Speech and Language Processing, Vol. 16, No. 1, 47-56, Jan. 2008.
%
% slight modifications by Eric Hunter to export some of the data.
% 'vf_segm' are the sections marked with vocal fry in 4ms windows at 2 ms
% steps. then you can use sum(vf_segm) to get the total segments.  and the
% length of the file.  output: frame_step_vf*sum(vf_segm) = total fry time
%
% addapted as a function 2018-Aug
%   filename - wave file to analyze
%   power_th - for normalized audio, this is about -20



% file name; power threshold (dB)
% filename = 'gdFANY04mCO1hD03.wav'; power_th = -20;
% filename = 'gdFANY04mCL1ae2s.wav'; power_th = -20;
% filename = 'gdFANY04mCL1ae3s.wav'; power_th = -20;
% filename = 'gdFANY04mCO3hA06.wav'; power_th = -20;
% filename = 'gdFANY04mCO3hE07.wav'; power_th = -20;
% filename = 'gdFANY04mCO3uu3s.wav'; power_th = -20;
% filename = 'gdFANY04mCL1ae3s110.wav'; power_th = -20;
% filename = 'gdFANY04mCL1ae3s441.wav'; power_th = -20;
% filename = 'gdFANY04mCL1ae3s000.wav'; power_th = -20;


% parameters for vocal fry detection
frame_size_vf = 0.004; % frame size for very short-term power estimation = 4 ms
frame_step_vf = 0.002; % frame step for very short-term power estimation = 2 ms
power_rise_th = 6; % minimum power rise = 6 dB
power_fall_th = 6; % minimum power fall = 6 dB
frame_size = 0.032; % frame size for short-term analysis = 32 ms
frame_step = 0.01; % 10 ms
ifp_th = 0.5; % maximum inter-frame periodicity = 0.5
ips_th = 0.5; % minimum inter-pulse similarity = 0.5
ipi_th = 0.12; % maximum inter-pulse interval = 120 ms
ips_frame_size = 0.005; % frame size for IPS estimation = 5 ms


% read waveform
[s,fs] = audioread(filename);
%fs = 16000; 
    

%% normalize and check sampling

if fs ~= 16000 % sampling frequency    
    s=resample(s,16000,fs);
%     soundsc(s,16000)
    fs=16000;
    fprintf('.fryFS')
%     disp('resampled audio')
end

s=0.95*s/max(abs(s));
% soundsc(s,16000)
%%

% parameter unit conversion for processing
frame_size = floor(fs * frame_size); % samples
frame_step = floor(fs * frame_step); % samples
frame_size_vf = floor(fs * frame_size_vf); % samples
frame_step_vf = floor(fs * frame_step_vf); % samples
ipi_th = floor( ipi_th * fs / frame_step_vf); % vst frames
ips_frame_size = floor(fs * ips_frame_size); % samples

% band-pass filtering
b = fir1(512,[100/fs 1500/fs],'bandpass');
s1 = filtfilt(b,1,s);

% (short-term) power estimation
power = [];
for n = 1:frame_step:length(s1)
    frame_initial_sample = floor(n - frame_size/2);

    for i = 1:frame_size
        if ( ((frame_initial_sample+i)>0) && ((frame_initial_sample+i)<=length(s1)) )
            frame(i) = s1(frame_initial_sample+i); % dc-removed orginal signal
        else
            frame(i) = 0;
        end
    end
    
    frame_power = 0; 
    for i=1:frame_size
        frame_power = frame_power + frame(i)*frame(i);
    end
    if (frame_power > 0)    
        frame_power = 10*log10(frame_power);
    else
        frame_power = -100;
    end
    
    frame_index = (n-1)/frame_step+1;
    power(frame_index) = frame_power;
end
nframes = frame_index;
time = [0:nframes-1]/fs*frame_step; % seconds


% very-short-term power estimation
vspower = [];
for n = 1:frame_step_vf:length(s1)
    frame_initial_sample = n - frame_size_vf/2;
    
    for i = 1:frame_size_vf
        if ( ((frame_initial_sample+i)>0) && ((frame_initial_sample+i)<=length(s1)) )
            frame(i) = s1(frame_initial_sample+i); % dc-removed orginal signal
        else
            frame(i) = 0;
        end
    end
    
    frame_power_vf = 0; 
    for i=1:frame_size_vf
        frame_power_vf = frame_power_vf + frame(i)*frame(i);
    end
    if (frame_power_vf > 0)
        frame_power_vf = 10*log10(frame_power_vf);
    else
        frame_power_vf = -100;
    end
    
    frame_index = (n-1)/frame_step_vf+1;
    vspower(frame_index) = frame_power_vf;
end
nframes_vf = frame_index;
time_vf = [0:nframes_vf-1]/fs*frame_step_vf; % seconds

% intra-frame periodicity estimation
ifp = [];
for n = 1:frame_step:length(s1)
    frame_initial_sample = floor(n - frame_size/2);

    % framing
    frame = [];
    for i = 1:frame_size
        if ( ((frame_initial_sample+i)>0) && ((frame_initial_sample+i)<=length(s1)) )
            %frame(i) = w(i)*s1(frame_initial_sample+i);
            frame(i) = s1(frame_initial_sample+i);
        else
            frame(i) = 0;
        end
    end
    
    % autocorrelation function
    for tau=0:frame_size-1
        acf(tau+1) = 0;
        for i=0:frame_size-1-tau
            acf(tau+1) = acf(tau+1) + frame(i+1)*frame(i+1+tau);
        end
    end
    if (acf(1)>0)
        for i = 2:frame_size
            acf(i) = acf(i) / acf(1);
        end
    end
    acf(1)=1;

    % detect peak in acf
    maxpeak = 0;
    maxpeaklag = 1;
    for i = 20:round(0.015*fs) % 15 ms
        if (acf(i)>maxpeak) && (acf(i)>acf(i-1)) && (acf(i)>acf(i+1))
            maxpeaklag = i;
            maxpeak = acf(i) * frame_size/(frame_size-maxpeaklag);
        end
    end
    
    % check sub-harmonics
    if maxpeaklag > 1
        k=2;
        j = k*maxpeaklag;
        while j < round(0.015*fs)
            maxpeak_sub = acf(j);
            for m=-3:3
                if acf(j+m)>acf(j)
                    maxpeak_sub = acf(j+m);
                    maxpeaklag_sub = j+m;
                end
            end
            maxpeak_sub = maxpeak_sub * frame_size/(frame_size-maxpeaklag_sub);
            if (maxpeak_sub < maxpeak)
                maxpeak = maxpeak_sub;
            end
            k=k+1;
            j = k*maxpeaklag;
        end
    end

    frame_index = (n-1)/frame_step+1;
    ifp(frame_index) = maxpeak;
end

% vocal fry pulse candidates
peak = [];
peak(nframes_vf) = 0;
for n = 6:nframes_vf-6
    if (vspower(n)-vspower(n-5) > power_rise_th) && (vspower(n)-vspower(n+5) > power_fall_th)
        % power peak detected
        peakpos = n;
        maxpower = vspower(n);
        for i=-5:5 % fine search
            if (vspower(n+i)>maxpower)
                vspower(n+i) = maxpower;
                peakpos = n+i;
            end
        end
        peak(peakpos) = 1;
        
        % corresponding frame index on short-term features
        n2 = round(n/(frame_step/frame_step_vf));
        if (n2 <= 0) 
            n2 = 1;
        end
       
        if ifp(n2) >= ifp_th % periodic
            peak(peakpos) = -0.2; % removed periodic candidate
        end
        if vspower(n) < power_th % low power
            peak(peakpos) = -0.1;
        end
            
        n = peakpos+5;
    end
end

% remove isolated peak candidates
for n = 6:nframes_vf-6
    if peak(n) > 0 % vf peak candidates
        flag_found = 0;
        for m = -ipi_th:ipi_th
            if (n+m > 0) && (n+m < nframes_vf) && (m~=0)
                if peak(n+m)>0
                    flag_found = 1;
                    break;
                end
            end
        end
        if flag_found == 0 % remove isolated candidates
            peak(n) = -0.3;
        end
    end
end

% estimate inter-pulse similarity for vocal fry pulse candidates
for n = 6:nframes_vf-6
    if peak(n) > 0 % vf peak candidates
        % estimate inter-pulse similarity
        ips = 0;
        for m = -ipi_th:-4 % closest peak is about 8 ms (4*2ms) away from the test peak
            if (n+m>0)&&(peak(n+m)>0)
                % compute ips with all peaks within the maximum allowed
                % inter-pulse interval
                for p = -0.005*fs/2:0.005*fs/2 
                    % compute cross-correlation 5 ms around the peak 
                    ccor = 0;
                    den1 = 0;
                    den2 = 0;
                    for q = -ips_frame_size/2:ips_frame_size/2
                        nsamp = n*frame_step_vf;
                        msamp = m*frame_step_vf;
                        ccor = ccor + s1(nsamp+p+q)*s1(nsamp+msamp+q);
                        den1 = den1 + s1(nsamp+p+q)*s1(nsamp+p+q);
                        den2 = den2 + s1(nsamp+msamp+q)*s1(nsamp+msamp+q);
                    end
                    if (den1>0) && (den2>0)
                        ccor = ccor / (sqrt(den1*den2));
                    else
                        ccor = 0;
                    end
                    if (ips<ccor) 
                        ips = ccor; % keep the largest ips
                    end
                end
                peak(n) = ips;
                if (peak(n+m)<ips)||(peak(n+m)==1)
                    peak(n+m) = ips;
                end
            end
        end
    end
end

% merge vocal fry segments
vf_segm = [];
vf_segm(nframes_vf) = 0;
for n = 6:nframes_vf-6
    if peak(n) > 0 % vf peaks
        for m = -ipi_th:-5
            if (n+m > 0) && (peak(n+m)>0)
                % merge
                for j = n+m+1:n-1
                    if peak(j)<=0
                        vf_segm(j) = 1;
                    end
                end
                vf_segm(n) = 1;
                vf_segm(n+m) = 1;
            end
        end
    end
end
                
% pre-roll (5ms)
for n = 6:nframes_vf-6
    if (vf_segm(n) > 0) && (vf_segm(n-1)<=0) 
         for j = 1:round(0.005*fs/frame_step_vf)
             vf_segm(n-j)=1;
         end
    end
end
% after-roll (10ms)
for n = nframes_vf-6:-1:6
    if (vf_segm(n) > 0) && (vf_segm(n+1)<=0) 
         for j = 1:round(0.01*fs/frame_step_vf)
             vf_segm(n+j)=1;
         end
    end
end

% vocal fry segment time
m = 1;
starttime = [];
endtime = [];
for n = 6:nframes_vf-6
    if (vf_segm(n) > 0) && (vf_segm(n-1)<=0) 
        % start time
        starttime(m) = n*frame_step_vf/fs;
    end
    if (vf_segm(n) > 0) && (vf_segm(n+1)<=0) 
        % end time
        endtime(m) = n*frame_step_vf/fs;
        m = m+1;
    end
end
vf_segm=logical(vf_segm);

%% plot results
if fignum>0
    figure(fignum),clf
    subplot(4,1,1)
    time_wave = [0:length(s1)-1]/fs;
    plot(time_wave,s1);
    axis([0 time_wave(length(s1)) min(s1) max(s1)]);
    title(filename);
    
    subplot(4,1,2)
    plot(time_vf,vspower,'b',[time_vf(1) time_vf(end)],[power_th power_th],'r');
    axis([0 time_vf(nframes_vf) max(min(vspower),-80) max(vspower)]);
    title('Very-short term power, power threshold (red)');
    ylabel('power (dB)');
    
    subplot(4,1,3)
    plot(time,ifp, [time_vf(1) time_vf(end)],[ifp_th ifp_th],'r');
    hold off;
    axis([0 time(nframes) 0 1]);
    title('Inter-frame periodicity (IFP), ifp threshold (red)');
    
    subplot(4,1,4)
    plot(time_vf,peak, time_vf,vf_segm,'g',[time_vf(1) time_vf(end)],[ips_th ips_th],'r');
    axis([0 time_vf(nframes_vf) -0.3 1.1]);
    title('Inter-pulse similarity (IPS)');
    xlabel('time (sec)')
    text(0.1, 0.4, ['fry time (sec): ' num2str(sum(vf_segm)*frame_step_vf/fs) ' sec'])
    
    saveas(gcf,[ '_fig_fryanalysis_' filename '.fig'])
    saveas(gcf,[ '_tif_fryanalysis_' filename '.tif'],'tiff')
    
end


% disp(' ')
% disp('--vocal fry analysis--')
% 
% disp(['filename, fry time (sec)'])
% disp([filename ', ' num2str(sum(vf_segm)*frame_step_vf/fs)])


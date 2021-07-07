function [outdata,sigaS,t1S,sig_dB,t2,t2Si,dB_vals,f0_So_aud,f0_fo_shp,f0_fo_prtf,f0_fo_aud,...
   g_duration] = a_9d4_par_task_function(...
   g_iter,g_cnt,g_filename,g_titlename,...
   s_gen,s_chnnum,s_gennum,s_avqi,s_gain_voltA,s_code,s_files,...
   s_RMS_window,s_foS_step,s_foS_Analysis,s_foS_NCand,s_foS_Accuracy,s_foS_SilenceThrsh,...
   s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost,...
   s_pertS_MaxPeriodFact,s_pertS_MaxAmpFact,...
   u_trimYN,u_concatYN,u_contour2xlsYN,u_fftenvelope,u_figson,u_showfindFofig,...
   c_LTAS_lower,c_LTAS_upper,c_autofo,c_foS_low,c_foS_up,c_fryanaly)

outdata=[]; sigaS=[]; t1S=[]; sig_dB=[]; t2=[]; t2Si=[]; dB_vals=[]; f0_So_aud=[]; f0_fo_shp=[]; f0_fo_prtf=[]; f0_fo_aud=[]; g_duration=[];

tic
%% 1_______________ Load new recording File _____________
%     j_iter = 1;
g_fname=g_filename{g_iter}; % get the filename of the current file.
fprintf('\n----- load %d of %d: %s -----',g_iter, g_cnt,g_fname)
fprintf('\n%dof%d',g_iter, g_cnt)
tmp = audioinfo([s_files '\' g_fname]);
g_Fs = tmp.SampleRate;
g_duration = tmp.Duration;
g_chn = tmp.NumChannels;

[siga,t1,g_Fs,g_fname]=LoadNewWave4(g_fname,g_Fs,[0 g_duration],s_files);
if g_chn>1    % two channel or one channel
    if s_chnnum==1  % channel 1, usually the acc on the collar
        siga=siga(:,1);
    elseif s_chnnum==2  % channel 2, usually the mic on the collar
        siga=siga(:,2);
    else   % maximum of the 2
        temprms1=rms(siga(:,1));
        temprms2=rms(siga(:,2));
        if temprms1>temprms2
            siga=siga(:,1);
        else
            siga=siga(:,2);
        end
    end
end
clear sig_RMS sig_dB sig_fftE t2 gain_volt f0_time3 f0_value3 f0_value1 f0_time1 s0_value1


%% set file-specific pitch ranges
if s_gen == 3
   if g_fname(s_gennum) == 'f' || g_fname(s_gennum) == 'F'
      s_g0 = 1;
      c_foS_low=150;    %lower bound for Fo (70-90 for males, 150 for females)
      c_foS_up=800;   %upper bound for Fo (350 for males, 800 for females)
   elseif g_fname(s_gennum) == 'm' || g_fname(s_gennum) == 'M'
      s_g0 = 2;
      c_foS_low=65;    %lower bound for Fo (70-90 for males, 150 for females)
      c_foS_up=350;   %upper bound for Fo (350 for males, 800 for females)
   end
end

c = fir1(34,[c_foS_low c_foS_up]/(g_Fs/2),'high');  % filter
sig = filtfilt(c,1,siga);
clear c
% sig = sig-mean(sig);        % zero


%% 1_______________ Get RMS and dB _____________
tmptoc=toc; fprintf(',dB')
[~,sig_dB,~,t2,~]=GetRMSandDBandSpecEnergy2(...
   siga,s_RMS_window,s_foS_step,g_Fs,s_gain_voltA);
t2=t2+t1(1);
if u_figson==1||u_figson==3
   figure(1),clf,subplot(3,1,1),cla
   plot(t2,sig_dB,'k',t1,siga/max(siga)*max(sig_dB),'r')
   axis([floor(t2(1)) t2(end) min(sig_dB) max(sig_dB)])
   legend('dB','sig','Location','southwest')
   axis([floor(t2(1)) t2(end) -max(sig_dB) max(sig_dB)])
   tmp=strrep([g_fname ', ' g_titlename ],'_', ' ');
   title([tmp ': [' num2str(g_iter) ' of ' num2str(g_cnt) ']' ],'Fontsize',8)
   drawnow
end

if u_trimYN>0
   try dB_clust=findclusters(sig_dB'); catch, dB_clust=findclusters(sig_dB(sig_dB>1)');end
   dB_thresh=min([dB_clust.low.stats.mean+1.5*dB_clust.low.stats.std ...
      dB_clust.high.stats.mean-1.5*dB_clust.high.stats.std]); % find the minimum
   g_fryfill=(10^((dB_clust.low.stats.mean-2*dB_clust.low.stats.std)/20))/s_gain_voltA;
else
   dB_thresh=sig_dB(1);
   g_fryfill=floor(mean([1 sig_dB(1)/2-1]));
   dB_clust.high.stats=basicstats(sig_dB);
   dB_clust.low.stats=basicstats(sig_dB);
end

clear cluster1 cluster1stats cluster2 cluster2stats obj idx tmp
fprintf('=%.0f',toc-tmptoc)
%_2_______________ Get RMS and dB _____________


%% 1_______________ Get Vocal Fry Areas _____________
tmptoc=toc; fprintf(',fry')
sigaMs=[];sigaFs=[];t1F2=[];t1M2=[];sigaM=[];t1M=[];sigaF=[];t1F=[]; %setup in case there is no fry

try
   [time_vf,vf_segm]=fun_vocalfrydetection_ejh(g_fname,-25,0);
   %        stp1=round(median(diff(t2))/(median(diff(time_vf))));
catch
   time_vf=0; vf_segm=true;
end
if sum(vf_segm)>2
   c_tmpFry=single(interp1(time_vf,double(vf_segm),t1,'PCHIP')');
   sigaFi=(c_tmpFry>0);                                  % logical of what is Fry
   sigaF=double(c_tmpFry).*siga; sigaFs=sigaF(sigaFi);    % shortened Fry only signal
   t1F=t1(sigaFi);t1F2=(1:length(t1F))/g_Fs;             % time to go with shortened signal
   sigaM=double(~c_tmpFry).*siga; sigaMs=sigaM(~sigaFi);
   t1M=t1(~sigaFi);t1M2=(1:length(t1M))/g_Fs;             % time to go with shortened signal
else
   c_tmpFry=siga.*0;
   sigaFi=logical(c_tmpFry);
end

tmp=basicstats(abs(siga));
tmp=g_fryfill*rand(length(siga),1)*tmp.Q1/2-tmp.Q1/4;tmp=tmp/6;

%    if u_contour2xlsYN==1
%       % write out the title and then the data in a single line
%       g_fid_out = fopen([g_titlename '_contour.csv'],'a'); %open file and create fid
%       fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'ModFry','time_vf');
%       dlmwrite([g_titlename '_contour.csv'], time_vf, '-append');
%       fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'ModFry','vf_segm');
%       dlmwrite([g_titlename '_contour.csv'], vf_segm, '-append');
%       fclose(g_fid_out); %close file
%    end

if c_fryanaly == 1 % replace modal voice as signal
   siga=double(~c_tmpFry).*siga+double(c_tmpFry).*tmp;
   sig=double(~c_tmpFry).*sig+double(c_tmpFry).*tmp;
end
if c_fryanaly == 2 % replace fry voice as signal
   siga=double(c_tmpFry).*siga+double(~c_tmpFry).*tmp;
   sig=double(c_tmpFry).*sig+double(~c_tmpFry).*tmp;
end

if sum(sigaFi)>0
   [tmptmM, ~] = praat_pitchGen3(0.95*sigaMs/max(abs(sigaMs)),...
      g_Fs,s_code,'praat_pitchGen3.psc',s_foS_Analysis,...
      s_foS_step/1000,30,700,s_foS_NCand,s_foS_Accuracy,s_foS_SilenceThrsh,...
      s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost);
   c_fryperc=100*(sum(sigaFi)/g_Fs)/(length(tmptmM)*s_foS_step/1000);
else
   c_fryperc=0;
end

fprintf('=%.0f',toc-tmptoc)
%_2_______________ Get Vocal Fry Areas  _____________


%% 1_______________ Get RMS and dB _____________
tmptoc=toc; fprintf(',re.dB')
[~,sig_dB,~,t2,~]=GetRMSandDBandSpecEnergy2(...
   siga,s_RMS_window,s_foS_step,g_Fs,s_gain_voltA);
t2=t2+t1(1); sig_dB(sig_dB<0)=min(sig_dB(sig_dB>0));
[sig_RMSm,sig_dBm,sig_fftEm,t2m,~]=GetRMSandDBandSpecEnergy2(...
   sigaM,s_RMS_window,s_foS_step,g_Fs,s_gain_voltA);
[sig_RMSf,sig_dBf,sig_fftEf,t2f,~]=GetRMSandDBandSpecEnergy2(...
   sigaF,s_RMS_window,s_foS_step,g_Fs,s_gain_voltA);

if u_figson==1||u_figson==3
   figure(1),clf,subplot(3,1,1),cla
   plot(t2,sig_dB,'k',t1,siga/max(siga)*max(sig_dB)/2+max(sig_dB)/3,'r',t2f,sig_dBf,'m')
   legend('dB','sig','fry','Location','southwest')
   axis([floor(t2(1)) t2(end) 0 max(sig_dB)])
   tmp=strrep([g_fname ', ' g_titlename ],'_', ' ');
   title([tmp ': [' num2str(g_iter) ' of ' num2str(g_cnt) ']' ],'Fontsize',9)
   figure(3), clf, subplot(6,1,1) % fo plot
   plot(t2,sig_dB,'k',t1,siga/max(siga)*max(sig_dB)/2+max(sig_dB)/3,'r',t2f,sig_dBf,'m')
   title([tmp ': [' num2str(g_iter) ' of ' num2str(g_cnt) ']' ],'Fontsize',9)
   axis([floor(t2(1)) t2(end) 0 max(sig_dB)]), drawnow
   clear tmp
end

%    if u_contour2xlsYN==1
%       % write out the title and then the data in a single line
%       g_fid_out = fopen([g_titlename '_contour.csv'],'a'); %open file and create fid
%       while g_fid_out<0
%          fprintf('...waiting...')
%          g_fid_out = fopen([g_titlename '_contour.csv'],'a'); %open file and create fid
%          fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b')
%       end
%       fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'sig_dB','t2'); %write out title
%       dlmwrite([g_titlename '_contour.csv'], t2, '-append'); % write out data
%       fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'sig_dB','dB');
%       dlmwrite([g_titlename '_contour.csv'], sig_dB, '-append');
%       fclose(g_fid_out); %close file
%    end
%    if u_contour2xlsYN==1
%       % write out the title and then the data in a single line
%       g_fid_out = fopen([g_titlename '_contour.csv'],'a'); %open file and create fid
%       fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'dBmodfry','t2'); %write out title
%       dlmwrite([g_titlename '_contour.csv'], t2, '-append'); % write out data
%       fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'dBmodfry','dB');
%       dlmwrite([g_titlename '_contour.csv'], sig_dB, '-append');
%       fclose(g_fid_out); %close file
%    end
fprintf('=%.0f',toc-tmptoc)
%_2_______________ Get RMS and dB _____________


%% 1_______________ Activity Level _________________
fprintf(',act')
ASLdata= ASL(siga,g_Fs); % use default parameters
c_ActiveSpeechLevel=ASLdata.activeSpeechLevel;
c_activityFactor=ASLdata.activityFactor;
clear ASLdata
%    fprintf('=%.0f',toc-tmptoc)
%_2 _______________ Activity Level


%% 1_______________ Estimate Pitch Range _____________
if c_autofo ==1
   tmptoc=toc; fprintf(',FoFnd')
   c_foS_lowi=c_foS_low;
   c_foS_upi=c_foS_up;
   [freqRef0,freqRef1,freqRef2,c_swapFo4Aud] = getFOrange2(siga,t1,sig_dB,t2,g_Fs,s_code,u_figson,u_showfindFofig,...
      c_foS_lowi,c_foS_upi,'praat_pitchGen3.psc',s_foS_Analysis,s_foS_step,s_foS_NCand,s_foS_Accuracy,...
      s_foS_SilenceThrsh,s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost);
   if c_swapFo4Aud==0
      c_foS_lowi=double(floor(freqRef1(1)/5)*5);
      c_foS_upi=double(ceil(freqRef1(2)/5)*5);
   else
      c_foS_lowi=double(floor(freqRef2(1)/5)*5);
      c_foS_upi=double(ceil(freqRef2(2)/5)*5);
      fprintf('**prob**FoFnd')
   end
   fprintf('=%.0f',toc-tmptoc)
   if (u_figson==1||u_figson==3)&&u_showfindFofig>0
      figure(u_showfindFofig)
      saveas(gcf,[s_files '\_tif' g_titlename '_' g_fname '_fofind.tif'],'tiff')
   end
else
      c_foS_lowi=c_foS_low;   
   c_foS_lowi=c_foS_low;
   c_foS_upi=c_foS_up;
end
% 2_______________ Estimate Pitch Range _____________


%% 1_______________ Estimate Pitch using PRAAT_____________
tmptoc=toc; fprintf(',PRT')
try
   [f0_tm_prt, f0_fo_prt] = praat_pitchGen3(0.95*siga/max(abs(siga)),g_Fs,s_code,'praat_pitchGen3.psc',...
      s_foS_Analysis,s_foS_step/1000,c_foS_lowi,c_foS_upi,s_foS_NCand,s_foS_Accuracy,s_foS_SilenceThrsh,...
      s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost);
catch
   fprintf(',rePRT')
   [f0_tm_prt, f0_fo_prt] = praat_pitchGen3(0.95*siga/max(abs(siga)),g_Fs,s_code,'praat_pitchGen3.psc',...
      'cc',s_foS_step/1000,c_foS_lowi,c_foS_upi,s_foS_NCand,s_foS_Accuracy,s_foS_SilenceThrsh,...
      s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost);
end
[f0_fo_prtf] = MakeSameLength(f0_fo_prt,f0_tm_prt,t2); % full size praat fo
if u_figson==1||u_figson==3
   figure(3),subplot(3,1,2), plot(f0_tm_prt,f0_fo_prt,'.r'),hold on
end
fprintf('=%.0f',toc-tmptoc)
% 2_______________ Estimate Pitch using PRAAT_____________


%% 1_______________ Estimate Pitch using SHRP _____________
tmptoc=toc; fprintf(',SHR')
tmp=max([c_foS_upi 3*c_foS_lowi]);
[temptime,f0_fo_shp,f0_SHR_shp,~]=shrp(0.95*siga/max(abs(siga)),...
   g_Fs,[c_foS_lowi tmp],6*s_foS_step,s_foS_step,0.1,tmp,20,0);
f0_tm_shp=temptime/1000; clear temptime
[f0_fo_shp] = MakeSameLength(f0_fo_shp,f0_tm_shp,t2);
[f0_SHR_shp] = MakeSameLength(f0_SHR_shp,f0_tm_shp,t2); clear f0_tm_shp
fprintf('=%.0f',toc-tmptoc)
if u_figson==1||u_figson==3
   figure(3),subplot(3,1,2), plot(t2,f0_fo_shp,'.g'),
end
clear tmp
%  2_______________ Estimate Pitch using SHRP _____________


%% 1_______________ Estimate Pitch using AudSwipe _____________
tmptoc=toc; fprintf(',Adswp')
[f0_fo_aud,f0_tm_aud,f0_So_aud] = audswipep(siga,g_Fs,[c_foS_lowi c_foS_upi],s_foS_step/1000,1/48,.1,0.5,0.05);
f0_So_aud=100*f0_So_aud; f0_So_aud(f0_So_aud<0)=0;
if u_figson==1||u_figson==3
   figure(3),subplot(3,1,2), plot(t2,f0_fo_aud,'.b'), hold off
   ylabel('f_0(hz)'),legend('praat', 'shrp','audswipe''','Location','best'),%legend('boxoff')
   xlabel('time (sec)')
   axis([t2(1) t2(end) c_foS_lowi c_foS_upi ])
   subplot(6,1,2)
   plot(t2,f0_So_aud,'.r'), legend('P_s audswipe''','Location','best')
   ylabel('P_s'), %legend('boxoff')
   axis([t2(1) t2(end) 0 80])
   drawnow
end
fprintf('=%.0f',toc-tmptoc)
% 2_______________ Estimate Pitch using AudSwipe _____________


%% 1_______________ Speech Segmentation & Detection _____________
fprintf(',Vc')
if c_fryanaly~=2
   if u_trimYN>0
      % dB clustering sort
      %          try dB_clust=findclusters(sig_dB'); catch dB_clust=findclusters(sig_dB(sig_dB>1)');end
      dB_vals = floor(min(sig_dB)):1:ceil(max(sig_dB));
      sig_dB_noise_dist = pdf(dB_clust.low.gcluster,dB_vals);
      sig_dB_speech_dist = pdf(dB_clust.high.gcluster,dB_vals);
      tmp=3:length(sig_dB_noise_dist)-3;
      [~, tmpi] = min(abs(sig_dB_noise_dist(tmp)./sig_dB_speech_dist(tmp)-1)); %find intersection of the two pdf
      tmp=dB_clust.high.stats.mean-2.7*dB_clust.high.stats.std;
      dB_thresh=floor(mean(dB_vals(tmpi))); % find the minimum
      % minimum of the overlap point compared to 2.7*std below mean
      if u_trimYN==4
          tmp=dB_thresh;
          try
              if dB_clust.low.stats.n/dB_clust.high.stats.n<0.1||dB_clust.high.stats.std<4 %only one cluster
                  dB_thresh=dB_vals(1);
              end
          catch
              if dB_clust.high.stats.std<4 % very narrow
                  dB_thresh=dB_vals(1);
              end
          end
          
          if u_figson==1||u_figson==3
              tmpdB=hist(sig_dB,dB_vals);tmpdB=2*tmpdB/sum(tmpdB);
              figure(3),subplot(3,1,3),
              plot(dB_vals,sig_dB_noise_dist,'ro-',dB_vals,sig_dB_speech_dist,'bs-',dB_vals,tmpdB,'xg-',...
                  [dB_thresh dB_thresh],[min([sig_dB_noise_dist sig_dB_speech_dist]) max(sig_dB_speech_dist)],'m-',...
                  [tmp tmp],[min([sig_dB_noise_dist sig_dB_speech_dist]) max(sig_dB_speech_dist)],'m--'), xlabel('dB')
              legend({'nois','spch'},'FontSize',6,'Location','northwest'),legend('boxoff')
              clear tmpdB
          end
          
      else
          
          % https://courses.lumenlearning.com/boundless-statistics/chapter/the-normal-curve/
          
          if u_figson==1||u_figson==3
              tmpdB=hist(sig_dB,dB_vals);tmpdB=2*tmpdB/sum(tmpdB);
              figure(3),subplot(3,1,3),
              plot(dB_vals,sig_dB_noise_dist,'ro-',dB_vals,sig_dB_speech_dist,'bs-',dB_vals,tmpdB,'xg-',...
                  [dB_thresh dB_thresh],[min([sig_dB_noise_dist sig_dB_speech_dist]) max(sig_dB_speech_dist)],'m-',...
                  [dB_vals(tmpi) dB_vals(tmpi)],[min([sig_dB_noise_dist sig_dB_speech_dist]) max(sig_dB_speech_dist)],'m--',...
                  [tmp tmp],[min([sig_dB_noise_dist sig_dB_speech_dist]) max(sig_dB_speech_dist)],'m--'), xlabel('dB')
              legend({'nois','spch'},'FontSize',6,'Location','northwest'),legend('boxoff')
              clear tmpdB
          end
      end
      
      tmpdB=sig_dB' > dB_thresh;
      tmpdBi=find(tmpdB>0);
      t2Si=logical(0.*sig_dB');
      t2Si(tmpdBi(1):tmpdBi(end))=1;
      tmp03F=single(interp1(t2,double(t2Si),t1,'PCHIP')');
      tmp01=logical(t2*0);tmp02=logical(t2*0);
      t2Si2=t2Si; % where dB is over the threshold
      
      tmpdB=sig_dB' > dB_thresh; % find dB above threshold
      tmpP0=f0_fo_aud>0;         % find fo from audswipe above 0
      tmpS0=f0_So_aud>.05;       % find so from audswiope above low pitch strength
      tmpF0=f0_fo_prtf>0;        % find fo from praat above 0
      
      if u_trimYN==2
         % get some thresholds,
         tmpI=tmpP0+tmpF0+tmpdB;
         tmpI=tmpI>=2;  % find where there is agreement between two of the three
         tmp=tmpI(1:end-1)+tmpI(2:end); tmp=tmp>1;
         tmp=tmp(1:end-1)+tmp(2:end);tmp=tmp>0;tmp=[0 tmp' 0]';
         if tmpI(1)+tmpI(2)==2, tmp(1)=1; end  % if the first 2 instances were voice
         if tmpI(end-1)+tmpI(end)==2, tmp(end)=1; end  % if the first 2 instances were voice
         tmpI=logical(tmp);
         f0_fo_shp=f0_fo_shp.*tmpI;  % SHRP with same voicing as PRAAT and Audswipe
         f0_fo_shp(f0_fo_shp==0)=NaN.*f0_fo_shp(f0_fo_shp==0);
         
         % do a running average both forward and backwards to smooth odd
         % ones and breaks out.
         windowSize = floor(.0005*g_Fs);
         tmp01=filter(ones(1,windowSize)/windowSize,1,tmpI);
         tmp02=filter(ones(1,windowSize)/windowSize,1,flipud(tmpI));
         tmp02=flipud(tmp02);
         t2Si2=(tmp01+tmp02).*tmpdB; %combine them
         t2Si=(t2Si2>.3);    % get threshold, this is a logical on what to keep
         tmp03F=single(interp1(t2,double(t2Si),t1,'PCHIP')');
      end
      
      if u_trimYN==3 %Fo only extraction
         tmpI=tmpP0+tmpF0;
         tmpI=tmpI>=1;  % find where there is agreement between the 2
         
         tmp=tmpI(1:end-1)+tmpI(2:end); tmp=tmp>1;
         tmp=tmp(1:end-1)+tmp(2:end);tmp=tmp>0;tmp=[0 tmp' 0]';
         if tmpI(1)+tmpI(2)==2, tmp(1)=1; end  % if the first 2 instances were voice
         if tmpI(end-1)+tmpI(end)==2, tmp(end)=1; end  % if the first 2 instances were voice
         tmpI=logical(tmp);
         
         t2Si2=logical(tmpI.*tmpdB);
         t2Si=t2Si2;tmp01=t2Si;tmp02=t2Si;
         tmp03F=single(interp1(t2,double(t2Si),t1,'PCHIP')');
      end
      
      if u_trimYN==4  %% get the end poinst and do it again
         % get some thresholds,
         tmpI=(tmpP0|tmpF0)&tmpdB;
         tmp=tmpI(1:end-1)+tmpI(2:end); tmp=tmp>1;
         tmp=tmp(1:end-1)+tmp(2:end);tmp=tmp>0;tmp=[0 tmp' 0]';
         if tmpI(1)+tmpI(2)==2, tmp(1)=1; end  % if the first 2 instances were voice
         if tmpI(end-1)+tmpI(end)==2, tmp(end)=1; end  % if the first 2 instances were voice
         tmpI=logical(tmp);
         tmp=find(tmpI>0);tmpI=0.*tmpI; tmpI(tmp(1):tmp(end))=1;
         f0_fo_shp=f0_fo_shp.*tmpI;  % SHRP with same voicing as PRAAT and Audswipe
         f0_fo_shp(f0_fo_shp==0)=NaN.*f0_fo_shp(f0_fo_shp==0);
         t2Si=logical(tmpI);
         tmp03F=single(interp1(t2,double(t2Si),t1,'PCHIP')');
      end
   else %keep whole signal
      t2Si=logical(1.*sig_dB');t2Si2=t2Si;tmpdB=t2Si;
      tmp01=logical(t2*0)';tmp02=logical(t2*0)';
      tmp03F=logical(siga*0+1);
      dB_stats=basicstats(sig_dB);
      dB_thresh=sig_dB(sig_dB>min([dB_stats.Q1 dB_stats.Q3]));
      clear dB_stats
   end
   if u_figson==1||u_figson==3
      figure(4)
      subplot(3,1,2),
      plot(t2,t2Si2,'yx',t2,tmp01,'bo',t2,tmp02,'ro',t2,t2Si,'g.',t2,tmpdB,'k.')
      title('looking for speech activity')
      subplot(3,1,1),
      plot(t1,tmp03F.*siga,'.',t1,tmp03F*max(siga),'g')
      xlabel('time')
      title('speech estimated area')
   end
   sigaSi=(tmp03F>0);                             % logical of what is voicing
   sigaS=double(tmp03F).*siga; sigaS=sigaS(sigaSi);       % shortened speech only signal
   t1S=t1(sigaSi);t1S2=(1:length(t1S))/g_Fs;        % time to go with shortened signal
else
   if sum(logical(c_tmpFry))>g_Fs/10
      sigaS=siga(logical(c_tmpFry));
      if length(sigaS)<2^(131)
         tmp=2^(13)-length(sigaS);
         sigaS=[0; 0; sigaS; zeros(tmp-2,1)];
         t1S=t1(logical(c_tmpFry));
         t1S=[t1S t1S(end)+(1:(tmp))/g_Fs];
         t1S2=(1:length(sigaS))/g_Fs;
      else
         t1S=t1(logical(c_tmpFry));t1S2=(1:length(sigaS))/g_Fs;
      end
   else
      sigaS=NaN;
      t1S=0;
   end
   sigaSi=sigaFi;
   s=resample(double(sigaFi),100,g_Fs);
   t2Si=abs(s)>0; t2Si=t2Si'; clear s
   f0_fo_prtf=NaN*t2Si';
   tmp=basicstats(sig_dB(~t2Si));
   dB_thresh=tmp.mean;
end
%      soundsc(sigaS,Fs)
clear tmp tmp01 tmp02 tmp03F tmpdB tmpdBi tmpF0 tmpfoM tmpi tmpI tmpP0 tmpS0 tmptmM
%  2_______________ Speech Segmentation & Detection _____________


%% 1_______________ pull some basic info for outputs
tmptoc=toc; fprintf(',FoStts')
fig=2;if u_figson==1||u_figson==3, figure(fig), clf, end
% all fo extractions in one
[FoAllhz_stats,FoAllst_stats] = freq_stats([f0_fo_prtf(t2Si)' f0_fo_aud(t2Si)' f0_fo_shp(t2Si)']',...
   [t2(t2Si) t2(t2Si) t2(t2Si)]',220);
tmp=[f0_fo_prtf(t2Si)' f0_fo_aud(t2Si)' f0_fo_shp(t2Si)']';
tmp=tmp(~isnan(tmp)); tmpt=t2(t2Si);
tmpt=[tmpt tmpt tmpt]'; tmpt=tmpt(~isnan(tmp)); p = polyfit(tmpt,tmp,1);
FoAllhz_stats.slope=p;
if u_fftenvelope>0
   tmpt=t2(t2Si); tmpt=[tmpt tmpt+t2(end) tmpt+2*t2(end)]';
   tmpt=tmpt(~isnan(tmp));  tmpt2=[t2 t2+t2(end) t2+2*t2(end)];
   yy = pchip([0 tmpt'],[mean(tmp) tmp'],[t2 t2+t2(end) t2+2*t2(end)]); tmp2=find(tmp>0);
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,2), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(yy(tmp2(1):tmp2(end)),tmpt2(tmp2(1):tmp2(end)),4,u_figson,'fo all');
   if ~isnan(linfitpnts),if u_figson==1||u_figson==3,figure(fig),subplot(3,3,2),axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]),end,end
   FoAllhz_stats.fft1stpeaksHzSlope=linfitpnts;
   FoAllhz_stats.fft1stpeaksHz=frq_out(:,2)';
   FoAllhz_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   FoAllhz_stats.fft1stpeaksHzSlope=NaN;
   FoAllhz_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   FoAllhz_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end

tmp=f0_So_aud(t2Si);tmp=tmp(~isnan(tmp)); % AudSwipe Ps
tmpt=t2(t2Si)'; tmpt=tmpt(~isnan(tmp)); p = polyfit(tmpt,tmp,1);
SoAud_stats= basicstats(round(tmp)); SoAud_stats.slope=p;
if u_fftenvelope>0
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,4), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(f0_So_aud',t2,4,u_figson,'So aud');
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,4), axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]), end
   SoAud_stats.fft1stpeaksHzSlope=linfitpnts;  % the slope of the first 4 peaks in dB envelope spectrum
   SoAud_stats.fft1stpeaksHz=frq_out(:,2)';
   SoAud_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   SoAud_stats.fft1stpeaksHzSlope=NaN;
   SoAud_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   SoAud_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end

dB_stats= basicstats(sig_dB);
p = polyfit(t2,sig_dB,1); dB_stats.slope=p;
if u_fftenvelope>0
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,5), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(sig_dB,t2,4,u_figson,'dB');
   if ~isnan(linfitpnts),if u_figson==1||u_figson==3,figure(fig),subplot(3,3,5),axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]),end,end
   dB_stats.fft1stpeaksHzSlope=linfitpnts;  % the slope of the first 4 peaks in dB envelope spectrum
   dB_stats.fft1stpeaksHz=frq_out(:,2)';
   dB_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   dB_stats.fft1stpeaksHzSlope=NaN;
   dB_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   dB_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end

tmp=sig_dB(t2Si)';tmp=tmp(~isnan(tmp)); % dB where in the upper distribution
tmpt=t2(t2Si)'; tmpt=tmpt(~isnan(tmp)); p = polyfit(tmpt,tmp,1);
dBv_stats= basicstats(tmp);  dBv_stats.slope=p;
if u_fftenvelope>0
   yy = pchip([0 tmpt'],[mean(tmp) tmp'],t2); tmp2=find(t2Si>0);
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,6), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(yy(tmp2(1):tmp2(end)),t2(tmp2(1):tmp2(end)),4,u_figson,'dBv');
   if ~isnan(linfitpnts),if u_figson==1||u_figson==3,figure(fig),subplot(3,3,6),axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]),end,end
   dBv_stats.fft1stpeaksHzSlope=linfitpnts;  % the slope of the first 4 peaks in dB envelope spectrum
   dBv_stats.fft1stpeaksHz=frq_out(:,2)';
   dBv_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   dBv_stats.fft1stpeaksHzSlope=NaN;
   dBv_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   dBv_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end

[FoSHRhz_stats,FoSHRst_stats] = freq_stats(f0_fo_shp(t2Si),t2(t2Si)',220);
if u_fftenvelope>0
   tmp=~isnan(f0_fo_shp);tmp=logical(single(t2Si).*single(tmp));
   yy = pchip([0 t2(tmp)],[mean(f0_fo_shp(tmp)) f0_fo_shp(tmp)'],t2); tmp2=find(tmp>0);
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,7), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(yy(tmp2(1):tmp2(end)),t2(tmp2(1):tmp2(end)),4,u_figson,'fo shp');
   if ~isnan(linfitpnts),if u_figson==1||u_figson==3,figure(fig),subplot(3,3,7),axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]),end,end
   FoSHRhz_stats.fft1stpeaksHzSlope=linfitpnts;
   FoSHRhz_stats.fft1stpeaksHz=frq_out(:,2)';
   FoSHRhz_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   FoSHRhz_stats.fft1stpeaksHzSlope=NaN;
   FoSHRhz_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   FoSHRhz_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end


[FoPRThz_stats,FoPRTst_stats] = freq_stats(f0_fo_prtf(t2Si),t2(t2Si)',220);
tmp=~isnan(f0_fo_prtf);tmp=logical(single(t2Si).*single(tmp));
yy = pchip([0 t2(tmp)],[mean(f0_fo_prtf(tmp)) f0_fo_prtf(tmp)'],t2); tmp2=find(tmp>0);
if u_fftenvelope>0
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,8), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(yy(tmp2(1):tmp2(end)),t2(tmp2(1):tmp2(end)),4,u_figson,'fo prtf');
   if ~isnan(linfitpnts),if u_figson==1||u_figson==3,figure(fig),subplot(3,3,8),axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]),end,end
   FoPRThz_stats.fft1stpeaksHzSlope=linfitpnts;
   FoPRThz_stats.fft1stpeaksHz=frq_out(:,2)';
   FoPRThz_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   FoPRThz_stats.fft1stpeaksHzSlope=NaN;
   FoPRThz_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   FoPRThz_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end

[FoAUDhz_stats,FoAUDst_stats] = freq_stats(f0_fo_aud(t2Si),t2(t2Si)',220);
if u_fftenvelope>0
   tmp=~isnan(f0_fo_aud);tmp=logical(single(t2Si).*single(tmp));
   yy = pchip([0 t2(tmp)],[mean(f0_fo_aud(tmp)) f0_fo_aud(tmp)'],t2); tmp2=find(tmp>0);
   if u_figson==1||u_figson==3, figure(fig), subplot(3,3,9), end
   [frq_out, linfitpnts,~]=FFT_peaks_slopes(yy(tmp2(1):tmp2(end)),t2(tmp2(1):tmp2(end)),4,u_figson,'fo aud');
   if ~isnan(linfitpnts),if u_figson==1||u_figson==3,figure(fig),subplot(3,3,9),axis([0 3*frq_out(4,2) min([3*frq_out(4,1) -1]) 0]),end,end
   FoAUDhz_stats.fft1stpeaksHzSlope=linfitpnts;
   FoAUDhz_stats.fft1stpeaksHz=frq_out(:,2)';
   FoAUDhz_stats.fft1stpeaksHzA=frq_out(:,1)';
else
   FoAUDhz_stats.fft1stpeaksHzSlope=NaN;
   FoAUDhz_stats.fft1stpeaksHz=[NaN NaN NaN NaN];
   FoAUDhz_stats.fft1stpeaksHzA=[NaN NaN NaN NaN];
end

if u_fftenvelope>0
if u_figson==1||u_figson==3
   figure(fig)
   saveas(gcf,[s_files '\_tif' g_titlename '_' g_fname '_arrayspect.tif'],'tiff')
   figure(3)
   saveas(gcf,[s_files '\_tif' g_titlename '_' g_fname '_foextract.tif'],'tiff')
end
end

% duration of the sample
tmp=find(t2Si>0);
durOspch1=s_foS_step/1000+t2(tmp(end))-t2(tmp(1)); %time between last and first point per speech detection (time then includes pauses)
if u_trimYN>0, durOspch2=s_foS_step*sum(t2Si)/1000; else, durOspch2=NaN; end % amount of speech identified (pauses removed)
tmp=[f0_fo_prtf(t2Si)' f0_fo_aud(t2Si)' f0_fo_shp(t2Si)'];
durOvoice=s_foS_step*sum(~isnan(tmp))/3/1000; % voicing time, using all 3 extraction techniques and comparing mutual values
durOvoicePr=s_foS_step*length(f0_fo_aud(~isnan(f0_fo_aud)))/1000; % voicing time per Praat only estimated
perOvoice1= 100*durOvoice/durOspch2;   %voicing percentage compared to judged amount of speech (voicing time of the concatinated sections)
perOvoice2= 100*durOvoice/durOspch1;   %voicing percentage compared to the beginning and end of the speech start (trimmed but with pauses)

clear tmp tmpt tmp2 output linfitpnts yy fig
% 2_______________ pull some basic info for outputs


%% 1_______________replot the all pitch extraction methods
%     plot(single(t2Si).*f0_value1');
if u_figson==1||u_figson==3
   %       tmp=220*2.^([FoAllst_stats.mean-FoAllst_stats.std FoAllst_stats.mean+FoAllst_stats.std FoAllst_stats.mean]/12);
   figure(1),subplot(3,1,2),
   plot(t2,f0_fo_prtf,'.m',t2(t2Si),f0_fo_prtf(t2Si),'xr',t2,f0_fo_aud,'.g',...
      t2(t2Si),f0_fo_aud(t2Si),'xb',t2(t2Si),f0_fo_shp(t2Si),'.c',...
      [t2(1) t2(end)],[FoAllhz_stats.median FoAllhz_stats.median],'k--',...
      [t2(1) t2(end)],[FoAllhz_stats.Q1 FoAllhz_stats.Q1],'b:',...
      [t2(1) t2(end)],[FoAllhz_stats.Q3 FoAllhz_stats.Q3],'b:',...
      [t2(1) t2(end)],[FoAllhz_stats.stmean FoAllhz_stats.stmean],'r--',...
      [t2(1) t2(end)],[FoAllhz_stats.stmeanmstd FoAllhz_stats.stmeanmstd],'m:',...
      [t2(1) t2(end)],[FoAllhz_stats.stmeanpstd FoAllhz_stats.stmeanpstd],'m:',...
      [t2(1) t2(end)],[FoAllhz_stats.stmeanm2std FoAllhz_stats.stmeanm2std],'m:',...
      [t2(1) t2(end)],[FoAllhz_stats.stmeanp2std FoAllhz_stats.stmeanp2std],'m:',...
      [t2(1) t2(end)],[FoAllhz_stats.stmeanm3std FoAllhz_stats.stmeanm3std],'m:',...
      [t2(1) t2(end)],[FoAllhz_stats.stmeanp3std FoAllhz_stats.stmeanp3std],'m:')
   axis([floor(t2(1)) t2(end) min([FoAllhz_stats.stmeanm3std 0.95*c_foS_lowi]) max([1.05*c_foS_upi FoAllhz_stats.stmeanp3std])])
   legend1 = legend('prt','prt_{cln}','aud','aud_{cln}','SHP','Location','NorthWest');
   set(legend1,'Location','Best','FontSize',6); legend('boxoff'),ylabel('F_0 (Hz)')
   tmp=(round(10*[FoAllhz_stats.mean FoAllhz_stats.stmean])/10); xlabel('time (sec)')
   title(['f_0_{(med)}=' num2str(round(10*FoAllhz_stats.median)/10) 'hz; '...
      'f_0_{(IQR)}=' num2str(round(10*FoAllhz_stats.IQR)/10) 'hz; '...
      '\muf_0=' num2str(tmp(1)) 'hz; ' '\muf_0=' num2str(tmp(2)) 'hz(st); '],'Fontsize',8)
   
   figure(1),subplot(3,1,1)
   tmp=dB_clust.high.stats;%if tmp.IQR==0,tmp.IQR=.05*tmp.max; tmp.median=tmp.max;end
   rngtmp=siga/max(abs(siga));rngtmp=0.97*rngtmp.*tmp.mean+20;%fo_st2_up=fo_st2_up+0.97*tmp.median;
   tmp3=(1:5:length(t1))';
   tmp4=sig_dBf'.*0+0.90*min(rngtmp);tmp4(sig_dBf>0)=1.03*max(rngtmp);
   plot(t2,sig_dB,'m',t1(tmp3),rngtmp(tmp3),'k',...
      t2,sig_dB.*(t2Si'),'.r',t1(tmp3),rngtmp(tmp3).*sigaSi(tmp3),'c:',t2f,tmp4,'b--')
   axis([floor(t2(1)) t2(end) max([-30 min(rngtmp)]) 1.03*max(rngtmp)]), ylabel('dB')
   legend(['dB_{\sigma}=' num2str(round(10*dBv_stats.std)/10) 'dB' ]), legend('boxoff')
   tmp=strrep([g_fname ', ' g_titlename ],'_', ' '); title(tmp,'Fontsize',8), drawnow
   clear tmp3 tmp legend1
end

if u_contour2xlsYN==1
   g_fid_out = fopen([g_titlename '_contour.csv'],'a'); %open file and create fid
   while g_fid_out<0
      fprintf('...waiting...')
      g_fid_out = fopen([g_titlename '_contour.csv'],'a'); %open file and create fid
      fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b')
   end
   fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'sig_dB','t2'); %write out title
   dlmwrite([g_titlename '_contour.csv'], t2, '-append'); % write out data
   fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'sig_dB','dB');
   dlmwrite([g_titlename '_contour.csv'], sig_dB, '-append');
   fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'speech_indx','fo');
   dlmwrite([g_titlename '_contour.csv'], t2Si', '-append'); % write out data
   fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'fo_PRAAT','fo');
   dlmwrite([g_titlename '_contour.csv'], f0_fo_prtf', '-append'); % write out data
   fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'fo_AudSwp','fo');
   dlmwrite([g_titlename '_contour.csv'], f0_fo_aud', '-append'); % write out data
   fprintf(g_fid_out,'%s,\t%s,\t%s,', g_fname,'fo_SHRP','fo');
   dlmwrite([g_titlename '_contour.csv'], f0_fo_shp', '-append'); % write out data
   fclose(g_fid_out); %close file
end

fprintf('=%.0f',toc-tmptoc)
clear cluster1 gcluster1 cluster1stats cluster2 gcluster2 cluster2stats
clear tmp2 tmp4 tmpdBi obj idx tmp windowSize
% 2_______________ replot the two pitch extraction methods _____________


%% 1_______________ Jitter & Shimmer & Entropy measures at pitch intervals
% description of some of the commands is here: http://www.fon.hum.uva.nl/praat/manual/Voice_2__Jitter.html
tmptoc=toc; fprintf('\n ... %dof%d,pert',g_iter, g_cnt)
tmp=0.95*sigaS/max(abs(sigaS)); c_jitter1=[];
if length(tmp)>10
   try
      [c_jitter1] = praat_voiceGen3(tmp,g_Fs,s_code,'praat_voiceGen3.psc','cc',...
         0,c_foS_lowi,c_foS_upi,s_foS_NCand,'no',0.03,0.25,...
         s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost,s_pertS_MaxPeriodFact,s_pertS_MaxAmpFact);
      
      if ~isfield(c_jitter1, 'jitter_abs')
         tmptoc=toc; fprintf(',repert')
         tmp=0.95*[sigaS sigaS]/max(abs([sigaS sigaS]));
         [jitter2] = praat_voiceGen3(tmp,g_Fs,s_code,'praat_voiceGen3.psc','cc',...
            0,c_foS_lowi,c_foS_upi,s_foS_NCand,'no',0.03,0.25,...
            s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost,s_pertS_MaxPeriodFact,s_pertS_MaxAmpFact);
      end
   catch
      disp(' --- problem getting perturbation measures ---')
   end
end
clear tmp

if ~isfield(c_jitter1, 'jitter_abs')
   clear jitter1
   c_jitter1.jitter=NaN;c_jitter1.jitter_abs=NaN;
   c_jitter1.jitter_rap=NaN;c_jitter1.jitter_ddp=NaN;
   c_jitter1.shimmer=NaN;c_jitter1.shimmer_db=NaN;
   c_jitter1.shimmer_apq3=NaN;c_jitter1.shimmer_apq5=NaN;
   c_jitter1.shimmer_apq11=NaN;c_jitter1.shimmer_dda=NaN;
   c_jitter1.nhr=NaN;c_jitter1.hnr=NaN;
end
fprintf('=%.0f',toc-tmptoc)

tmptoc=toc; fprintf(',trpy')
% Get RPDE and DFA
tmp=resample(sigaS, g_Fs, 25000); tmp=tmp-mean(tmp); tmp=tmp/max(abs(tmp));
try %works only on 32 bit windows
   c_jitter1.rpde = rpde(tmp, 4, 50, 0.2, 1000);  % sig, d-Embedding dimension, tau-Embedding delay, eta-RPDE close returns radius, Tmax
   [dfa, ~, ~] = fastdfa(tmp, (50:20:200)'); % DFA scaling range 50-200
   c_jitter1.dfa= 1/(1+exp(-dfa)); %alpha_norm
catch
   tmptoc=toc; fprintf(',skpRDP&DFA')
   c_jitter1.rpde = NaN;
   c_jitter1.dfa=NaN;
end

try  % Get PPE
   c_jitter1.ppe = ppe2([f0_fo_prt ; f0_fo_aud(~isnan(f0_fo_aud))]);
catch
   c_jitter1.ppe =NaN;
   tmptoc=toc; fprintf(',skpPPE')
end
clear tmp

c_jitter1.length=length(sigaS)/g_Fs;
fprintf('=%.0f',toc-tmptoc)
% 2_______________ Jitter and Shimmer at pitch intervals


%% 1_______________ Mark S Analysis
tmptoc=toc; fprintf(',HFCC')
clear mark1 mark1E mark1parameters
HFCCparm.fsHFCC=16000;
HFCCparm.HPF_order=10;
HFCCparm.HPF_fc=5000;
HFCCparm.trimThreshold=25;
HFCCparm.trimTimeConstant=20e-3;
HFCCparm.calcASL=false;
HFCCtmp.x=(sigaS);
HFCCtmp.fs=g_Fs;
try
   [HFCCtmp,~] = processSentences(HFCCtmp,HFCCparm);
   c_HFCC.HFCC_SDS=HFCCtmp.ccMeasures.stDevSum;
   c_HFCC.DHFCC_SDS=HFCCtmp.DccMeasures.stDevSum;
   c_HFCC.DDHFCC_SDS=HFCCtmp.DDccMeasures.stDevSum;
catch
   c_HFCC.HFCC_SDS=NaN;
   c_HFCC.DHFCC_SDS=NaN;
   c_HFCC.DDHFCC_SDS=NaN;
end

clear HFCCtmp HFCCtmpE HFCCparm
fprintf('=%.0f',toc-tmptoc)
% 1_______________ Mark Analysis


%% 1_______________ LTAS and Alpha Ratio _____________
tmptoc=toc; fprintf(',LTAS')
if length(sigaS)>1
   try
      LTAS=ltas2(sigaS/max(abs(sigaS)),g_Fs,min([2^floor(log2(length(sigaS))) 2^13]),...
         'han',0,0,1,43,c_LTAS_lower,c_LTAS_upper,FoAllhz_stats.median);
   catch
      LTAS=ltas2([sigaS;sigaS;sigaS]/max(abs(sigaS)),g_Fs,min([2^floor(log2(length(sigaS))) 2^13]),...
         'han',0,0,1,43,c_LTAS_lower,c_LTAS_upper,FoAllhz_stats.median);
   end
   LTAS.sig=[];
   f1=LTAS.f;
   %       spec1dB=LTAS.norm0spectrum;
else
   LTAS.AlphaRatio=NaN;LTAS.dB1kTo3k=NaN;
   LTAS.dB1kTo3kN=NaN;LTAS.totdBfft=NaN;
   LTAS.dBspect_slope=NaN;
   LTAS.dBspect_slopeFo=NaN;
   LTAS.tilt = NaN;
end

if u_figson==1||u_figson==3
   tmp=1:length(LTAS.dBspectrumFlat);
   figure(1),subplot(3,1,3),
   plot(f1(tmp),LTAS.dBspectrumFlat,'m',...
      f1(tmp),LTAS.norm0spectrum(tmp)-LTAS.norm0spectrum(1),'g',...
      [  50 1000],-20*[1 1],'b.-.',...
      [1000 5000],-20*[1 1],'r.-.',...
      [1000 3150],-30*[1 1],'k.--',...
      [2800 3400],-40*[1 1],'b.-',...
      [c_LTAS_lower c_LTAS_upper],LTAS.dBspect_slope*[c_LTAS_lower c_LTAS_upper]+max(LTAS.dBspectrumFlat),'b--',...
      [FoAllhz_stats.median c_LTAS_upper],LTAS.dBspect_slopeFo*[FoAllhz_stats.median c_LTAS_upper]+max(LTAS.dBspectrumFlat),'g-.')
   ylabel('LTAS (dB)'),xlabel('Hz'),drawnow
   title(['\alpha_{ratio}=' num2str(round(LTAS.AlphaRatio)) ...
       'dB ; {\SigmadB_{1-3.125k}}=' num2str(round(10*LTAS.dB1kTo3k)/10) ...
       'dB,' num2str(100*round(1000*LTAS.dB1kTo3kN)/1000) ...
       '% ; {\SigmadB_{2.8-3.4k}}=' num2str(round(10*LTAS.dBsingersform)/10)...
       'dB,', num2str(100*round(1000*LTAS.dBsingersformN)/1000) '%' ],'Fontsize',8)
   %       text(100,-40,['\alpha_{ratio}=' num2str(round(LTAS.AlphaRatio)) ...
   %          'dB ; {\SigmadB_{n(1-3k)}}=' num2str(round(100*LTAS.dB1kTo3kN)/100) ...
   %          'dB ; {\SigmadB_{n(2.8-3.4k)}}=' num2str(round(100*LTAS.dBsingersformN)/100)] )
   %       text(100,-35,['slope=' num2str(round(10000*LTAS.dBspect_slope)/100) ...
   %          '; slope_{f0}=' num2str(round(10000*LTAS.dBspect_slopeFo)/100) ...
   %          '; tilt=' num2str(round(100*LTAS.tilt)/100)] )
   text(50,-20+3,'\alpha_{ratio}a')
   text(1000,-20+3,'\alpha_{ratio}b')
   text(1000,-30+3,'\Sigma(1-3.15kHz)')
   text(2800,-40+3,'SingForm')
   drawnow
   try
      axis([0 c_LTAS_upper ...
         min([-45 LTAS.dBspect_slope*[c_LTAS_lower c_LTAS_upper]+max(LTAS.dBspectrumFlat)])...
         max(LTAS.dBspectrumFlat)])
   catch
      text(c_LTAS_upper/2,0.5,'all NaN')
   end
   
end
%    clear tmp tmp1 tmp2
clear spec1 spec1dB f1 tmp
fprintf('=%.0f',toc-tmptoc)
%_2_______________ LTAS and Alpha Ratio_____________


%% CPPS / AVQI analysis
clear outavqiC outavqiO

[cppC,~,~] = CPP(0.99999999*sigaS/max(abs(sigaS)),g_Fs,[c_foS_lowi, c_foS_upi]);
[cppO,~,~] = CPP(0.99999999*siga/max(abs(siga)),g_Fs,[c_foS_lowi, c_foS_upi]);

if s_avqi == 1
   tmptoc=toc; fprintf(',CPPSv')
   % [outavqiC] = praat_cpps(sigaS/max(abs(sigaS)), g_Fs, g_fname, s_code,s_files);
   [outavqiC] = praat_cpps3(0.99999999*sigaS/max(abs(sigaS)), g_Fs,'praat_avqi3.praat', s_code, c_foS_lowi, c_foS_upi);
   fprintf('=%.0f',toc-tmptoc)
else
   outavqiC.cpps = NaN;
   outavqiC.hnr = NaN;
   outavqiC.shim = NaN;
   outavqiC.shdb = NaN;
   outavqiC.slope = NaN;
   outavqiC.tilt = NaN;
   outavqiC.avqi = NaN;
end
outavqiC.cpp=cppC;
if s_avqi == 1
   tmptoc=toc; fprintf(',CPPSall')
   [outavqiO] = praat_cpps3(0.99999999*siga/max(abs(siga)), g_Fs,'praat_avqi3.praat', s_code, c_foS_lowi, c_foS_upi);
   % [outavqiO] = praat_cpps(siga/max(abs(siga)), g_Fs, g_fname, s_code,s_files);
   fprintf('=%.0f',toc-tmptoc)
   if u_figson==1||u_figson==3
      figure(1),subplot(3,1,1),
      xlabel(['CPPS_v=' num2str(round(10*outavqiC.cpps)/10) ...
         'dB; CPPS_{all}=' num2str(round(10*outavqiO.cpps)/10) 'dB'])
   end
else
   outavqiO.cpps = NaN;
   outavqiO.hnr = NaN;
   outavqiO.shim = NaN;
   outavqiO.shdb = NaN;
   outavqiO.slope = NaN;
   outavqiO.tilt = NaN;
   outavqiO.avqi = NaN;
   fprintf(',skp.CPPS')
end
outavqiO.cpp=cppO;
clear cppC cppO

%%  Save plots
fprintf(',Figs')
if u_figson==1||u_figson==3
   figure(1)
   saveas(gcf,[s_files '\_fig' g_titlename '_' g_fname '.fig'])
   saveas(gcf,[s_files '\_tif' g_titlename '_' g_fname '.tif'],'tiff')
end
if u_concatYN==2
   tmp=[s_files '\' g_titlename '_reduced_' g_fname(1:end-4) '.wav'];
   audiowrite(tmp,concat1.sigS ,g_Fs)
end


%% write out summary
tmptoc=toc; fprintf(',Save')

% SETTING OUTPUT
%    u_trimYN                 TrimOptn
%    c_fryanaly               fryanaly
%    s_RMS_window             RMSwindw
%    s_foS_step               FoWinStp
%    s_foS_SilenceThrsh       ptSlncTh
%    s_foS_VoiceThrsh         ptVoicTh
%    s_foS_OctCost            prOctCst
%    s_foS_OctJumpCost        prOctJmp
%    s_foS_VoiceUnvoiceCost   prVcUnVc
%    s_pertS_MaxPeriodFact    prMxPerF
%    s_pertS_MaxAmpFact       prMxAmpF
%    c_foS_lowi               FoExtLwr
%    c_foS_upi                FoExtUp_
%    c_LTAS_lower             LTASlwr_
%    c_LTAS_upper             LTASuppr
% GENERAL DATA
%    durOspch1                DurBeg2End
%    durOspch2                DurSpchDtc
%    durOvoice                DurVoicALL
%    durOvoicePr              DurVoicPRT
%    perOvoice1               %VoicCnSpc
%    perOvoice2               %VoicAlSpc
%    g_duration               FileTime__
% GENERAL CALCULATED AND SPECTRAL DATA
%    c_activityFactor         ActyFact
%    c_fryperc                FryPerc
%    LTAS.totdBfft            spdB_fft
%    LTAS.totdBrms            spdB_rms
%    LTAS.AlphaRatio          AlphaRto
%    LTAS.dBspect_slope       dBspcSlp
%    LTAS.dBspect_slopeFo     dBspcSpF
%    LTAS.tilt                dBspcSpO
%    LTAS.dB1kTo3k            dB1kTo3k
%    LTAS.dB1kTo3kN           dB1kT3kN
%    LTAS.dBsingersform       dBSngFm
%    LTAS.dBsingersformN      dBSngFmN
%    c_HFCC.HFCC_SDS          HFCCSDS
%    c_HFCC.DHFCC_SDS         DHFCCSDS
%    c_HFCC.DDHFCC_SDS        DDHFCSDS
% CPP SPEECH ONLY
%    outavqiC.cpps            cpps_avC
%    outavqiC.hnr             hnr_avC
%    outavqiC.shim            shim_avC
%    outavqiC.shdb            shdb_avC
%    outavqiC.slope           slop_avC
%    outavqiC.tilt            tilt_avC
%    outavqiC.avqi            avqi_avC
% CPP WHOLE FILE
%    outavqiO.cpps            cpps_avO
%    outavqiO.hnr             hnr_avO
%    outavqiO.shim            shim_avO
%    outavqiO.shdb            shdb_avO
%    outavqiO.slope           slop_avO
%    outavqiO.tilt            tilt_avO
%    outavqiO.avqi            avqi_avO
% PURTERBATION MEASURES (ON SPEECH SEGMENT ONLY, NO PAUSES)
%    c_jitter1.jitter         jitter
%    c_jitter1.jitter_abs     jit_abs
%    c_jitter1.jitter_rap     jit_rap
%    c_jitter1.jitter_ddp     jit_ddp
%    c_jitter1.shimmer        shimmer
%    c_jitter1.shimmer_db     shim_db
%    c_jitter1.shimmer_apq3   shmapq3
%    c_jitter1.shimmer_apq5   shmapq5
%    c_jitter1.shimmer_apq11  shmapq11
%    c_jitter1.shimmer_dda    shim_dda
%    c_jitter1.nhr            nhr
%    c_jitter1.hnr            hnr
%    c_jitter1.rpde           rpde
%    c_jitter1.dfa            dfa
%    c_jitter1.ppe            ppe
%    c_jitter1.length         purt_dur
% SPEECH LEVEL OF THE CONNECTED SPEECH SEGMENTS (NO PAUSES)
%    SoAud_stats
%    FoAllhz_stats
%    FoAllst_stats
%    FoPRThz_stats
%    FoPRTst_stats
%    FoAUDhz_stats
%    FoAUDst_stats
%    FoSHRhz_stats
%    FoSHRst_stats

outdata=NaN;     % metrics
outdataD='----------';    % metrics names

% Parameters
outdata=[outdata,NaN,...
   u_trimYN,c_fryanaly,s_RMS_window,s_foS_step,...     %15
   s_foS_SilenceThrsh,s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,...
   s_foS_VoiceUnvoiceCost,s_pertS_MaxPeriodFact,s_pertS_MaxAmpFact,...
   c_swapFo4Aud,c_foS_lowi,c_foS_upi,c_LTAS_lower,c_LTAS_upper];
outdataD=[outdataD;'*SETTINGS*';...
   'trimOption'; 'fry_analy_'; 'RMS_windw_'; 'Fo_WinStp_';...
   'perSilncTh'; 'perVoicTh_'; 'perOctCst_'; 'perOctJmp_';...
   'perVcUnVcC'; 'perMxPerF_'; 'perMxAmpF_'; ...
   'Fo_ExtProb'; 'Fo_LwrBndr'; 'Fo_uprBndr'; 'LTASlwr_Bd'; 'LTASupprBd'];%slope accross utterance
outdataD=[outdataD;'----------'];
outdata=[outdata,NaN];

outdata=[outdata,NaN,...     %5
   durOspch1,...   % end voicing - beginning voicing
   durOspch2,...   % estimated length of all speaking parts
   durOvoice,...   % duration of voice only combining all 3 extraction types
   durOvoicePr,...   % duration of voice only from PRAAT raw
   perOvoice1,...   % voicing percentage compared to judged amount of speech (voicing time of the concatinated sections)
   perOvoice2,...   % %voicing percentage compared to the beginning and end of the speech start (trimmed but with pauses)
   g_duration];    % overall length of read file
outdataD=[outdataD;'*GEN_DATA_';...
   'DurBeg2End';...      % end voicing - beginning voicing
   'DurSpchDtc';...      % estimated length of all speaking parts
   'DurVoicALL';...      % duration of voice only (fo segments from all 3 extraction methods)
   'DurVoicPRT';...      % duration of voice only (Praat only voicing)
   '%VoicCnSpc';...      % Voice % compared to the concatenated speech
   '%VoicAlSpc';...      % Voice % from the length of speech, from first instance to last instance
   'FileTime__'];
outdataD=[outdataD;'----------'];
outdata=[outdata,NaN];

outdata=[outdata,NaN,...     %15
   c_activityFactor,...
   c_fryperc,...        % percent of recording which is fry
   LTAS.AlphaRatio,LTAS.dBspect_slope,LTAS.dBspect_slopeFo,LTAS.tilt,...
   LTAS.dB1kTo3k,...
   LTAS.dB1kTo3kN,...
   LTAS.dBsingersform,...
   LTAS.dBsingersformN,...
   c_HFCC.HFCC_SDS,...
   c_HFCC.DHFCC_SDS,...
   c_HFCC.DDHFCC_SDS];
outdataD=[outdataD;'*GEN_CALC_';...
   'ActvtyFact';...      % activity factor activityFactor [0 1]
   'FryPerc   ';...      % amount of the recording which is fry compared to the modal segments
   'AlphaRatio'; 'dBspecSlp_'; 'dBspecSlpF'; 'dBspcSpOct';... % alpha ratio, spec slope from 50hz, spec slope from Fo, slop using 3rd oct
   'dB1kTo3k  ';...      % energy from 1kHz to 3kHz normed to that range
   'dB1kTo3kN ';...      % percent energy from 1kHz to 3kHz compared to 50hz-10khz
   'dBSngFmnt ';...      % energy from 2.8kHz to 3.4kHz (singers formant area)
   'dBSngFmntN';...      % percent energy from 2.8kHz to 3.4kHz (singers formant area) to LTAS area
   'HFCCSDS   ';...      % HFCC HFCC standard deviation sum measure
   'DHFCCSDS  ';...      % delta HFCC features, standard deviation sum measure
   'DDHFCSDS  '];        % delta delta HFCC features, standard deviation sum measure
outdataD=[outdataD;'----------'];
outdata=[outdata,NaN];

outdata=[outdata,NaN,...     %7
   outavqiC.cpp,outavqiC.cpps,outavqiC.hnr,outavqiC.shim,outavqiC.shdb,...
   outavqiC.slope,outavqiC.tilt,outavqiC.avqi];
outdataD=[outdataD;'*CPPctSPCH';...
   'cpp__avCAT';'cpps_avCAT';'hnr_avCAT ';'shim_avCAT';'shdb_avCAT';...
   'slop_avCAT';'tilt_avCAT';'avqi_avCAT'];
outdataD=[outdataD;'----------'];
outdata=[outdata,NaN];

outdata=[outdata,NaN,...     %7
   outavqiO.cpp,outavqiO.cpps,outavqiO.hnr,outavqiO.shim,outavqiO.shdb,...
   outavqiO.slope,outavqiO.tilt,outavqiO.avqi];
outdataD=[outdataD;'*CPPallFIL';...
   'cpp__avALL';'cpps_avALL';'hnr_avALL ';'shim_avALL';'shdb_avALL';...
   'slop_avALL';'tilt_avALL';'avqi_avALL'];
outdataD=[outdataD;'----------'];
outdata=[outdata,NaN];

outdata=[outdata,NaN,...     %16
   c_jitter1.jitter,c_jitter1.jitter_abs,c_jitter1.jitter_rap,c_jitter1.jitter_ppq5,c_jitter1.jitter_ddp,...
   c_jitter1.shimmer,c_jitter1.shimmer_db,c_jitter1.shimmer_apq3,c_jitter1.shimmer_apq5,...
   c_jitter1.shimmer_apq11,c_jitter1.shimmer_dda,...
   c_jitter1.nhr,c_jitter1.hnr,...
   c_jitter1.rpde,c_jitter1.dfa,c_jitter1.ppe,...
   c_jitter1.length];

outdataD=[outdataD;'*CALC_PURT';...
   'jitter    ';'jitt_abs  ';'jitt_rap  ';'jit_ppq5  ';'jit_ddp   ';...
   'shimmer   ';'shim_db   ';'shim_apq3 ';'shim_apq5 ';...
   'shim_apq11';'shim_dda  ';...
   'nhr       ';'hnr       ';...
   'rpde      ';'dfa       ';'ppe       ';...
   'purt_dur  '];
outdataD=[outdataD;'----------'];
outdata=[outdata,NaN];

outdata=[outdata,NaN,...
   dB_stats.n,dB_stats.mean,dB_stats.std,dB_stats.var,...     %13
   dB_stats.mode,dB_stats.median,dB_stats.Q1,dB_stats.Q3,...
   dB_stats.IQR,dB_stats.skew,dB_stats.bowley_skew,dB_stats.kurt,...
   dB_stats.slope(1),dB_stats.fft1stpeaksHzSlope,dB_stats.fft1stpeaksHz(1:4),...
   dB_stats.fft1stpeaksHzA(1:4)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*dB_ALLraw';...
   'dB__N     '; 'dB__mean  '; 'dB__std   '; 'dB__var   ';...
   'dB__mode  '; 'dB__med   '; 'dB__Q1    '; 'dB__Q3    ';...
   'dB__IQR   '; 'dB__skew  '; 'dB__bowskw'; 'dB__kurt  ';...
   'dB__slop  '; 'dB__fftslp'; 'dB__ffthz1'; 'dB__ffthz2'; 'dB__ffthz3'; 'dB__ffthz4';...
   'dB__fftA1 '; 'dB__fftA2 '; 'dB__fftA3 '; 'dB__fftA4 '];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   dBv_stats.n,dBv_stats.mean,dBv_stats.std,dBv_stats.var,...     %13
   dBv_stats.mode,dBv_stats.median,dBv_stats.Q1,dBv_stats.Q3,...
   dBv_stats.IQR,dBv_stats.skew,dBv_stats.bowley_skew,dBv_stats.kurt,...
   dBv_stats.slope(1),dBv_stats.fft1stpeaksHzSlope,dBv_stats.fft1stpeaksHz(1:4),...
   dBv_stats.fft1stpeaksHzA(1:4)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*dBvALLspc';...
   'dBv_N     '; 'dBv_mean  '; 'dBv_std   '; 'dBv_var   ';...
   'dBv_mode  '; 'dBv_med   '; 'dBv_Q1    '; 'dBv_Q3    ';...
   'dBv_IQR   '; 'dBv_skew  '; 'dBv_bowskw'; 'dBv_kurt  ';...
   'dBv_slop  '; 'dBv_fftslp'; 'dBv_ffthz1'; 'dBv_ffthz2'; 'dBv_ffthz3'; 'dBv_ffthz4';...
   'dBv_Afft1 '; 'dBv_Afft2 '; 'dBv_Afft3 '; 'dBv_Afft4 '];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   SoAud_stats.n,SoAud_stats.mean,SoAud_stats.std,SoAud_stats.var,...     %13
   SoAud_stats.mode,SoAud_stats.median,SoAud_stats.Q1,SoAud_stats.Q3,...
   SoAud_stats.IQR,SoAud_stats.skew,SoAud_stats.bowley_skew,SoAud_stats.kurt,...
   SoAud_stats.slope(1),SoAud_stats.fft1stpeaksHzSlope,...
   SoAud_stats.fft1stpeaksHz(1:4),SoAud_stats.fft1stpeaksHzA(1:4)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Sov_PchSt';...
   'Sov_N     '; 'Sov_mean  '; 'Sov_std   '; 'Sov_var   ';...
   'Sov_mode  '; 'Sov_med   '; 'Sov_Q1    '; 'Sov_Q3    ';...
   'Sov_IQR   '; 'Sov_skew  '; 'Sov_bowskw'; 'Sov_kurt  ';...
   'Sov_slop  '; 'Sov_fftslp'; 'Sov_ffthz1'; 'Sov_ffthz2'; 'Sov_ffthz3'; 'Sov_ffthz4';...
   'Sov_Afft1 '; 'Sov_Afft2 '; 'Sov_Afft3 '; 'Sov_Afft4 '];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoAllhz_stats.n,FoAllhz_stats.mean,FoAllhz_stats.std,FoAllhz_stats.var,...     %13
   FoAllhz_stats.mode,FoAllhz_stats.median,FoAllhz_stats.Q1,FoAllhz_stats.Q3,...
   FoAllhz_stats.IQR,FoAllhz_stats.skew,FoAllhz_stats.bowley_skew,FoAllhz_stats.kurt,...
   FoAllhz_stats.slope(1),FoAllhz_stats.fft1stpeaksHzSlope,FoAllhz_stats.fft1stpeaksHz(1:4),...
   FoAllhz_stats.fft1stpeaksHzA(1:4),FoAllhz_stats.stmean,...
   FoAllhz_stats.stmeanpstd,FoAllhz_stats.stmeanmstd];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Fo_Hz_ALL';...
   'FoAl_N    '; 'FoAl_mean '; 'FoAl_std  '; 'FoAl_var  ';...
   'FoAl_mode '; 'FoAl_med  '; 'FoAl_Q1   '; 'FoAl_Q3   ';...
   'FoAl_IQR  '; 'FoAl_skew '; 'FoAlbowskw'; 'FoAl_kurt ';...
   'FoAl_slop '; 'FoAlfftslp'; 'FoAlffthz1'; 'FoAlffthz2'; 'FoAlffthz3'; 'FoAlffthz4';...
   'FoAl_Afft1'; 'FoAl_Afft2'; 'FoAl_Afft3'; 'FoAl_Afft4';...
   'FoAlmeanST'; 'FoAmnpsdST'; 'FoAmnmstST'];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoAllst_stats.n,FoAllst_stats.mean,FoAllst_stats.std,FoAllst_stats.var,...     %13
   FoAllst_stats.mode,FoAllst_stats.median,FoAllst_stats.Q1,FoAllst_stats.Q3,...
   FoAllst_stats.IQR,FoAllst_stats.skew,FoAllst_stats.bowley_skew,FoAllst_stats.kurt,...
   FoAllst_stats.slope(1)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Fo_ST_ALL';...
   'FoAl_N    '; 'FoAl_mean '; 'FoAl_std  '; 'FoAl_var  ';...
   'FoAl_mode '; 'FoAl_med  '; 'FoAl_Q1   '; 'FoAl_Q3   ';...
   'FoAl_IQR  '; 'FoAl_skew '; 'FoAlbowskw'; 'FoAl_kurt ';...
   'FoAl_slop '];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoPRThz_stats.n,FoPRThz_stats.mean,FoPRThz_stats.std,FoPRThz_stats.var,...     %13
   FoPRThz_stats.mode,FoPRThz_stats.median,FoPRThz_stats.Q1,FoPRThz_stats.Q3,...
   FoPRThz_stats.IQR,FoPRThz_stats.skew,FoPRThz_stats.bowley_skew,FoPRThz_stats.kurt,...
   FoPRThz_stats.slope(1),FoPRThz_stats.fft1stpeaksHzSlope,...
   FoPRThz_stats.fft1stpeaksHz(1:4),FoPRThz_stats.fft1stpeaksHzA(1:4),...
   FoPRThz_stats.stmean,...
   FoPRThz_stats.stmeanpstd,FoPRThz_stats.stmeanmstd];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Fov_Hz___';...
   'Fov_N     '; 'Fov_mean  '; 'Fov_std   '; 'Fov_var   ';...
   'Fov_mode  '; 'Fov_med   '; 'Fov_Q1    '; 'Fov_Q3    ';...
   'Fov_IQR   '; 'Fov_skew  '; 'Fov_bowskw'; 'Fov_kurt  ';...
   'Fov_slop  '; 'Fov_fftslp'; 'Fov_ffthz1'; 'Fov_ffthz2'; 'Fov_ffthz3'; 'Fov_ffthz4';...
   'Fov_Afft1 '; 'Fov_Afft2 '; 'Fov_Afft3 '; 'Fov_Afft4 ';...
   'Fov_meanST'; 'FovmnpsdST'; 'FovmnmsdST'];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoPRTst_stats.n,FoPRTst_stats.mean,FoPRTst_stats.std,FoPRTst_stats.var,...     %13
   FoPRTst_stats.mode,FoPRTst_stats.median,FoPRTst_stats.Q1,FoPRTst_stats.Q3,...
   FoPRTst_stats.IQR,FoPRTst_stats.skew,FoPRTst_stats.bowley_skew,FoPRTst_stats.kurt,...
   FoPRTst_stats.slope(1)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Fov_ST___';...
   'Fov_N     '; 'Fov_mean  '; 'Fov_std   '; 'Fov_var   ';...
   'Fov_mode  '; 'Fov_med   '; 'Fov_Q1    '; 'Fov_Q3    ';...
   'Fov_IQR   '; 'Fov_skew  '; 'Fov_bowskw'; 'Fov_kurt  ';...
   'Fov_slop  '];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoAUDhz_stats.n,FoAUDhz_stats.mean,FoAUDhz_stats.std,FoAUDhz_stats.var,...     %13
   FoAUDhz_stats.mode,FoAUDhz_stats.median,FoAUDhz_stats.Q1,FoAUDhz_stats.Q3,...
   FoAUDhz_stats.IQR,FoAUDhz_stats.skew,FoAUDhz_stats.bowley_skew,FoAUDhz_stats.kurt,...
   FoAUDhz_stats.slope(1),FoAUDhz_stats.fft1stpeaksHzSlope,...
   FoAUDhz_stats.fft1stpeaksHz(1:4),FoAUDhz_stats.fft1stpeaksHzA(1:4),...
   FoAUDhz_stats.stmean,...
   FoAUDhz_stats.stmeanpstd,FoAUDhz_stats.stmeanmstd];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Pov_Hz___';...
   'Pov_N     '; 'Pov_mean  '; 'Pov_std   '; 'Pov_var   ';...
   'Pov_mode  '; 'Pov_med   '; 'Pov_Q1    '; 'Pov_Q3    ';...
   'Pov_IQR   '; 'Pov_skew  '; 'Pov_bowskw'; 'Pov_kurt  ';...
   'Pov_slop  '; 'Pov_fftslp'; 'Pov_ffthz1'; 'Pov_ffthz2'; 'Pov_ffthz3'; 'Pov_ffthz4';...
   'Pov_Afft1 '; 'Pov_Afft2 '; 'Pov_Afft3 '; 'Pov_Afft4 ';...
   'Pov_meanST'; 'PovmnpsdST'; 'PovmnmsdST'];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoAUDst_stats.n,FoAUDst_stats.mean,FoAUDst_stats.std,FoAUDst_stats.var,...     %13
   FoAUDst_stats.mode,FoAUDst_stats.median,FoAUDst_stats.Q1,FoAUDst_stats.Q3,...
   FoAUDst_stats.IQR,FoAUDst_stats.skew,FoAUDst_stats.bowley_skew,FoAUDst_stats.kurt,...
   FoAUDst_stats.slope(1)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Pov_ST___';...
   'Pov_N     '; 'Pov_mean  '; 'Pov_std   '; 'Pov_var   ';...
   'Pov_mode  '; 'Pov_med   '; 'Pov_Q1    '; 'Pov_Q3    ';...
   'Pov_IQR   '; 'Pov_skew  '; 'Pov_bowskw'; 'Pov_kurt  ';...
   'Pov_slop  '];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoSHRhz_stats.n,FoSHRhz_stats.mean,FoSHRhz_stats.std,FoSHRhz_stats.var,...     %13
   FoSHRhz_stats.mode,FoSHRhz_stats.median,FoSHRhz_stats.Q1,FoSHRhz_stats.Q3,...
   FoSHRhz_stats.IQR,FoSHRhz_stats.skew,FoSHRhz_stats.bowley_skew,FoSHRhz_stats.kurt,...
   FoSHRhz_stats.slope(1),FoSHRhz_stats.fft1stpeaksHzSlope,...
   FoSHRhz_stats.fft1stpeaksHz(1:4),FoSHRhz_stats.fft1stpeaksHzA(1:4),...
   FoSHRhz_stats.stmean,...
   FoSHRhz_stats.stmeanpstd,FoSHRhz_stats.stmeanmstd];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Fo_SHR_Hz';...
   'FoS_N     '; 'FoS_mean  '; 'FoS_std   '; 'FoS_var   ';...
   'FoS_mode  '; 'FoS_med   '; 'FoS_Q1    '; 'FoS_Q3    ';...
   'FoS_IQR   '; 'FoS_skew  '; 'FoS_bowskw'; 'FoS_kurt  ';...
   'FoS_slop  '; 'FoS_fftslp'; 'FoS_ffthz1'; 'FoS_ffthz2'; 'FoS_ffthz3'; 'FoS_ffthz4';...
   'FoS_Afft1 '; 'FoS_Afft2 '; 'FoS_Afft3 '; 'FoS_Afft4 ';...
   'FoS_meanST'; 'FoSmnpsdST'; 'FoSmnmsdST'];%slope accross utterance
outdataD=[outdataD;'----------'];

outdata=[outdata,NaN,...
   FoSHRst_stats.n,FoSHRst_stats.mean,FoSHRst_stats.std,FoSHRst_stats.var,...     %13
   FoSHRst_stats.mode,FoSHRst_stats.median,FoSHRst_stats.Q1,FoSHRst_stats.Q3,...
   FoSHRst_stats.IQR,FoSHRst_stats.skew,FoSHRst_stats.bowley_skew,FoSHRst_stats.kurt,...
   FoSHRst_stats.slope(1)];
outdata=[outdata,NaN];
outdataD=[outdataD;'*Fo_SHR_ST';...
   'FoS_N     '; 'FoS_mean  '; 'FoS_std   '; 'FoS_var   ';...
   'FoS_mode  '; 'FoS_med   '; 'FoS_Q1    '; 'FoS_Q3    ';...
   'FoS_IQR   '; 'FoS_skew  '; 'FoS_bowskw'; 'FoS_kurt  ';...
   'FoS_slop  '];%slope accross utterance
outdataD=[outdataD;'----------'];

fprintf('=%.1f',toc)
%     keyboard

% outdataD

%% write out to success file
g_fid2_out = fopen([g_titlename '_trackcomplete.txt'],'a'); %open file and create fid
while g_fid2_out<0
   fprintf('...waiting...')
   g_fid2_out = fopen([g_titlename '_trackcomplete.txt'],'a'); %open file and create fid
end
fprintf(g_fid2_out,'%d\n',g_iter); %write out title
fclose(g_fid2_out); %close file




function [freqRef0,freqRef1,freqRef2,swapFo4Aud] = getFOrange2(siga,t1,sig_dB,t2,g_Fs,s_code,u_figson,fig,...
   c_foS_lowi,c_foS_upi,praatscrpt,s_foS_Analysis,s_foS_step,s_foS_NCand,s_foS_Accuracy,...
   s_foS_SilenceThrsh,s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% get estimate pitch range
% get dB for sample and segment into two sections
% try dB_clust=findclusters(sig_dB'); catch dB_clust=findclusters(sig_dB(sig_dB>1)');end
% dB_thresh=min([dB_clust.low.stats.mean+1.5*dB_clust.low.stats.std ...
%    dB_clust.high.stats.mean-1.5*dB_clust.high.stats.std]); % find the minimum
 dB_vals = floor(min(sig_dB)):1:ceil(max(sig_dB));
 
% get first pass fo
[f0_tm_prt, f0_hz_prt] = praat_pitchGen3(0.95*siga/max(abs(siga)),g_Fs,s_code,praatscrpt,...
   s_foS_Analysis,s_foS_step/1000,c_foS_lowi,c_foS_upi,s_foS_NCand,s_foS_Accuracy,s_foS_SilenceThrsh,...
   s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost);
[f0_hz_prt] = MakeSameLength(f0_hz_prt,f0_tm_prt,t2); clear f0_tm_prt
[f0_fo_aud,~,f0_So_aud] = audswipep(siga,g_Fs,[c_foS_lowi c_foS_upi],s_foS_step/1000,1/48,.1,0.5,0.05);
[temptime,f0_fo_shp,f0_SHR_shp,~]=shrp(0.95*siga/max(abs(siga)),...
   g_Fs,[c_foS_lowi c_foS_upi],6*s_foS_step,s_foS_step,0.1,c_foS_upi,20,0);
f0_tm_shp=temptime/1000;
[f0_fo_shp] = MakeSameLength(f0_fo_shp,f0_tm_shp,t2);
[f0_SHR_shp] = MakeSameLength(f0_SHR_shp,f0_tm_shp,t2); clear f0_tm_shp
 
%% speech segments should be where 2 of 4 exist: praat, audfo, ps>0.05, dB>dB thresh
dBi=sig_dB'>dB_vals(floor(length(dB_vals)/8)+1); % find where dB is above the threshold
Soi=f0_So_aud>0.05; % find where pitch strength is over 0.05
Poi=~isnan(f0_fo_aud); % find where there is voicing from audswipe
Foi=~isnan(f0_hz_prt); % find where there is voicing from praat
Spch1i=dBi+Soi+Poi+Foi; Spch1i=Spch1i>2;
%remove single instances
tmp=Spch1i(1:end-1)+Spch1i(2:end); tmp=tmp>1;
tmp=tmp(1:end-1)+tmp(2:end);tmp=tmp>0;tmp=[0 tmp' 0];
if Spch1i(1)+Spch1i(2)==2, tmp(1)=1; end  % if the first 2 instances were voice
if Spch1i(end-1)+Spch1i(end)==2, tmp(end)=1; end  % if the first 2 instances were voice
Spch1i=logical(tmp);
f0_fo_shpC=f0_fo_shp.*Spch1i';  % SHRP with same voicing as PRAAT and Audswipe
f0_fo_shpC(f0_fo_shpC==0)=NaN.*f0_fo_shpC(f0_fo_shpC==0);
So_thresh=0.1; So_i=f0_So_aud>So_thresh;  % pitch strength threshold for speech is about 0.1
SHR_thresh=0.05; SHR_i=f0_SHR_shp>SHR_thresh;

%% check if there is similarities...
% disp('fo stats')
tmp1=basicstats(f0_hz_prt(Foi)); 
% disp(['prt ' num2str(tmp1.mean) ', ' num2str(tmp1.std)])
tmp2=basicstats(f0_fo_aud(Poi)); 
% disp(['aud ' num2str(tmp2.mean) ', ' num2str(tmp2.std)])
tmp3=f0_fo_shpC(Spch1i); tmp3=tmp3(~isnan(tmp3));
tmp3=basicstats(tmp3); 
% disp(['shp ' num2str(tmp3.mean) ', ' num2str(tmp3.std)])

tmp4=std([tmp1.mean tmp2.mean tmp3.mean]);
tmp5=std([tmp1.std tmp2.std tmp3.std]);
swapFo4Aud=0;
if tmp4>30
%    disp('fo range problem')
   swapFo4Aud=1;
end

%% combine the extraction data for one output
t2_mult=[t2 t2 t2];
if swapFo4Aud==0
   fo_mult0=[ f0_hz_prt' f0_fo_aud' f0_fo_shpC'];  % combined Fo as they are
else
   fo_mult0=[ f0_fo_aud' f0_fo_aud' f0_fo_aud'];  % combined Fo as they are
end
SP_mult=[ Spch1i Spch1i Spch1i];
So_mult=[So_i' So_i' So_i'];
fo_mult1=fo_mult0.*SP_mult; fo_mult1(fo_mult1==0)=NaN; % combined Fo filtered with dB & combined voicing 
fo_mult2=fo_mult0.*SP_mult.*So_mult; fo_mult2(fo_mult2==0)=NaN; % same as fo_multi1 but with pitch strength threshold

% semitones, find the mode and then fit to it.
fo_x2_stats0= basicstats(fo_mult0(~isnan(fo_mult0))); % get non NaN fo's
fo_x2_stats1= basicstats(fo_mult1(~isnan(fo_mult1)));
fo_x2_stats2= basicstats(fo_mult2(~isnan(fo_mult2)));
fo_up0=fo_x2_stats0.mean+fo_x2_stats0.std;
fo_low0=fo_x2_stats0.mean-.5*fo_x2_stats0.std;
fo_up1=fo_x2_stats1.mean+fo_x2_stats1.std;
fo_low1=fo_x2_stats1.mean-.5*fo_x2_stats1.std;
fo_up2=fo_x2_stats2.mean+3*fo_x2_stats2.std;
fo_low2=fo_x2_stats2.mean-2*fo_x2_stats2.std;

fo_x2_st0 = 12*log2(fo_mult0./fo_x2_stats0.mean);
fo_x2_st1 = 12*log2(fo_mult1./fo_x2_stats1.mean);
fo_x2_st2 = 12*log2(fo_mult2./fo_x2_stats2.mean);

fo_x2_st0_stats= basicstats(fo_x2_st0(~isnan(fo_x2_st0)));
fo_st0_up=fo_x2_st0_stats.mean+3*fo_x2_st0_stats.std;
fo_st0_low=fo_x2_st0_stats.mean-2*fo_x2_st0_stats.std;

fo_x2_st1_stats= basicstats(fo_x2_st1(~isnan(fo_x2_st1)));
fo_st1_up=fo_x2_st1_stats.mean+3*fo_x2_st1_stats.std;
fo_st1_low=fo_x2_st1_stats.mean-2*fo_x2_st1_stats.std;

fo_x2_st2_stats= basicstats(fo_x2_st2(~isnan(fo_x2_st2)));
fo_st2_up=fo_x2_st2_stats.mean+3*fo_x2_st2_stats.std;
fo_st2_low=fo_x2_st2_stats.mean-2*fo_x2_st2_stats.std;

freqRef0= fo_x2_stats0.mean.*2.^([fo_st0_low fo_st0_up]/12);
freqRef1= fo_x2_stats1.mean.*2.^([fo_st1_low fo_st1_up]/12);
freqRef2= fo_x2_stats2.mean.*2.^([fo_st2_low fo_st2_up]/12);

%% figure plot
if (u_figson==1||u_figson==3)&&fig>0
   figure(fig),clf,subplot(10,1,1),
   plot(t2,sig_dB,'k',t1,siga/max(siga)*max(sig_dB),'r'),ylabel('dB')
   title(['starting range: ' num2str(c_foS_lowi) ' - ' num2str(c_foS_upi) 'hz' ])
   
   figure(fig),subplot(10,1,2), tmp=logical(Spch1i.*So_i');
   plot(t2,Spch1i,':b',t2,f0_So_aud,'oc',t2(Spch1i),f0_SHR_shp(Spch1i),'og',...
      t2(So_i),f0_So_aud(So_i),'.r',t2(tmp),f0_SHR_shp(tmp),'.m'),
   
   %    legend('spch_{area}','so','SHR','cln','FontSize',6),ylabel('segments')
   
   figure(fig),subplot(5,1,2),
   plot(t2,f0_fo_shp,'xg',t2,f0_hz_prt,'oc',t2,f0_fo_aud,'.m',t2_mult,fo_mult2,'.y',...
      [t2_mult(1) t2_mult(end)],[fo_up2 fo_up2],'--r',...
      [t2_mult(1) t2_mult(end)],[fo_low2 fo_low2],'--r'),
   ylabel('F_0 (Hz)'),%xlabel('time (sec)'),drawnow
%    legend({'fo_{shrp}','fo_{praat}','fo_{audswp}','combo'},'FontSize',6)
   
   figure(fig), subplot(5,3,7),hist(fo_mult0,40),xlabel('hz0'),ylabel('N')
   figure(fig), subplot(5,3,8),hist(fo_mult1,20),xlabel('hz1'),%ylabel('N')
   figure(fig), subplot(5,3,9),hist(fo_mult2,20),xlabel('hz2'),%ylabel('N')
   
   figure(fig),subplot(5,1,4),
   plot(t2_mult,fo_x2_st0,'og',t2_mult,fo_x2_st1,'.b',t2_mult,fo_x2_st2,'.r',...
      [t2_mult(1) t2_mult(end)],[fo_st2_low fo_st2_low],'--r',...
      [t2_mult(1) t2_mult(end)],[fo_st2_up fo_st2_up],'--r'),
   legend({'fo_{full}','fo_{spchi}','fo_{spchi+}'},'FontSize',8), legend BOXOFF
   ylabel('F_0 (ST)'),%xlabel('time (sec)'),drawnow
   
   figure(fig), subplot(5,3,13),hist(fo_x2_st0,40),xlabel('hz0'),ylabel('N')
   title(['low: ' num2str(round(freqRef0(1))) ', high ' num2str(round(freqRef0(2)))])
   figure(fig), subplot(5,3,14),hist(fo_x2_st1,20),xlabel('hz1')
   title(['low: ' num2str(round(freqRef1(1))) ', high ' num2str(round(freqRef1(2)))])
   figure(fig), subplot(5,3,15),hist(fo_x2_st2,20),xlabel('hz2')
   title(['low: ' num2str(round(freqRef2(1))) ', high ' num2str(round(freqRef2(2)))])
%    figure(fig),subplot(5,3,14), 
%    set(gca,'XTick',zeros(1,0),'XTickLabel',{});
%    set(gca,'YTick',zeros(1,0),'YTickLabel',{});
%    plot(10,10,'w.',0,0,'w.',0,0,'w.',0,0,'w.',0,0,'w.',0,0,'w.')
%    
%    set(gca,'XColor',[1 1 1],'XTick',zeros(1,0),'XTickLabel',{},'YColor',...
%       [1 1 1],'YTick',zeros(1,0),'YTickLabel',{},'ZColor',[1 1 1]);
%    % Create legend
%    legend1 = legend(['___CLEAN___'],...   %fname,...
%       ['low: ' num2str(round(freqRef2(1))) ', high ' num2str(round(freqRef2(2)))], ...
%       ['___ALL___'],...
%       ['low: ' num2str(round(freqRef0(1))) ', high: ' num2str(round(freqRef0(2)))],...
%       'Location','Best'); legend BOXOFF
%    set(legend1,'Interpreter','none','Location','north');
end


%%

end


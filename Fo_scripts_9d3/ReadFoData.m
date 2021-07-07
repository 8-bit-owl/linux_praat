

clear all;
clc;
% pTxt = ['\\ncvs-data2\NCVS_research\TobiasBirdSong\ZF1-24hours']; % add the path that has the data
% addpath(pTxt)
dB_overall = [];                                          % initialize the overall vectors
F0_overall = [];
time_overall = 0;
%% 1_______________ Load new recording File _____________
filename=FilenamesByExt('tmp.mat');
cnt=length(filename);

for n=1:cnt
  filename{n};
  
  load(filename{n})
  
  % adjust
  time=data4.time;
  dB=data4.dB;
  F0=data4.F0;
  fftE=data4.fftE;
  % clear data4
%   figure(10),plot(time,dB,'g.',time,F0,'r.')
  fig=1;
%   keyboard
  
  %%
  button='Redo';
  while button(1)=='R'
  figure(1),hist(dB,[min(dB):max(dB)])
  xlim([max([0 min(dB)]) max(dB)])
  title('get lower dB cutoff for voicing')
  tmp=ginput(1);
  dBlow=floor(tmp(1));
  indx1=find(dB>dBlow&F0<max(F0)&F0>min(F0));
  
  figure(1),hist((fftE),[min(fftE):max(fftE)])
  xlim([max([0 min(fftE)]) mean(fftE)+3*std(fftE)])
  title('get lower dB cutoff for voicing')
  tmp=ginput(1);
  dBElow=floor(tmp(1));
  indx2=find(fftE>dBElow&F0<max(F0)&F0>min(F0));

  tmp=F0*0;
  tmp(indx1)=1;tmp(indx2)=tmp(indx2)+1;
  indx=find(tmp==2);
%   
%   d_dB=[min(dB):max(dB)];
%   d_fftE=min(fftE): mean(fftE)+3*std(fftE);
%   dB_hist=hist(dB,d_dB);
%   dB_fftE=hist(fftE,d_fftE);
% 
%   plot(d_dB,dB_hist, d_fftE,dB_fftE),
%   
%   title('get lower dB cutoff for voicing')
%   tmp=ginput(1);
%   dBlow=floor(tmp(1));
%   indx=find(dB>dBlow&F0<max(F0)&F0>min(F0));
  
  figure(1),hist(F0(indx),[min(F0(indx)):5:max(F0(indx))])
  title('get lower F_0 cutoff for voicing')
  tmp=ginput(1);
  F0low=floor(tmp(1));F0low=max([F0low min(F0)]);
  
  title('get higher F_0 cutoff for voicing')
  tmp=ginput(1);
  F0high=floor(tmp(1));  
  hz=[80:6:1200];

  indx=find(fftE>dBElow&dB>dBlow&F0<F0high&F0>F0low);
  

%   indx=find(dB>40&F0<max(hz)&F0>min(hz));  
  FoStats=stats(F0(indx));
  
  % plots
%   figure(fig),fig=fig+1;
  figure(2)
  hist(F0(indx),hz)
  occ=hist(F0(indx),hz);
  title(filename{n})
  xlabel('F_0')
  ylabel('# occ')

  saveas(gcf,[ 'hist_' filename{n} '_bins' '.fig'])
  saveas(gcf,[ 'hist_' filename{n} '_bins' '.tif'],'tiff')  

  
  % hz=hz(2:end);occ=occ(2:end);
%   figure(fig),fig=fig+1;
  figure(3)
  plot(hz,occ)
  title(filename{n})  
  xlabel('F_0')
  ylabel('# occ')
  
  saveas(gcf,[ 'hist_' filename{n} '_line' '.fig'])
  saveas(gcf,[ 'hist_' filename{n} '_line' '.tif'],'tiff')  
  
  button = questdlg('Keep these results or Redo', ...
        'Exit Dialog','Keep','Redo','Keep');
  
  end
  
  FoStats1{n}=FoStats;
  FoHist.hz{n}=hz;
  FoHist.occ{n}=occ;
 
  plot(time(indx),F0(indx),'.')

  %%
  
end


%%
save('5min.mat', 'FoStats1', 'FoHist', 'FoHist','filename')

%%
keyboard

%% plot histogram
clear FoHist.occ_n
figure(1),hold on
for n=1:length(FoHist.hz)
  FoHist.occ_n(n,:)=FoHist.occ{n}/max(FoHist.occ{n});
  plot(Age(n),FoStats1{n}.mode,'b.')
end
% hold off
% plot(FoHist.hz{1},FoHist.occ_n)
% plot(FoHist.hz{1},mean(FoHist.occ_n))

% Ag


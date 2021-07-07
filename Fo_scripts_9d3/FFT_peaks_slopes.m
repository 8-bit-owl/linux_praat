function [output,linfit4pnts,curvfitall] = FFT_peaks_slopes(data,datatime,numpnts,fig,titlen)
%FFT_peaks_slopes take the fft of the incoming signal and find the slope of
%the spectrum as well as the first few peaks (and slope)
%   input
%     data - array of data of equal timesteps
%     datatime - time array going with data
%     varargin3 - plot graph and figure number

%%
Y=fft(data/max(data));
P2 = abs(Y/length(data)); 
P1 = P2(1:floor(length(data)/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
P1 = log10(P1);
P1=P1-max(P1(2:end));
[ymax,inx] = func_pickpeaks(P1(2:end), 5); 
if length(ymax)<numpnts
%    inx=1:numpnts; ymax=NaN*inx;
   output=ones(numpnts,2).*NaN;
   curvfitall=[NaN NaN NaN];
   linfit4pnts=NaN;
   if fig>0
      % figure(fig),
      plot(0,0,'rx')
      % title('Single-Sided Amplitude Spectrum')
      xlabel('f (Hz)')
      ylabel(['log(|' titlen '|)'])
   end
else
   inx=inx+1;
   f = 1/(mode(diff(datatime)))*(0:(length(data)/2))/length(data);
   output=[ymax(1:numpnts)' f(inx(1:numpnts))'];
   output=double(output);
   linfit4pnts = polyfit(f(inx(1:numpnts)),ymax(1:numpnts),1);
   linfit4pnts=linfit4pnts(1);
   
   if fig>0
      % figure(fig),
      plot(f(2:end),P1(2:end),f(inx),ymax,'ro')
      % title('Single-Sided Amplitude Spectrum')
      xlabel('f (Hz)')
      ylabel(['log(|' titlen '|)'])
      axis([f(1) f(end) 1.05*min(ymax) .05])
   end
   curvfitall = polyfit(f(inx),ymax,2);
end








function [ymax,inx] = func_pickpeaks(y, winlen)
% [ymax,inx] = func_pickpeaks(y, winlen)
% Input:  y - from wavread
%         winlen - size of window to use
% Output: ymax - peak values
%         inx - peak positions
% Notes:  
%
% Author: Yen-Liang Shue and Markus Iseli, Speech Processing and Auditory Perception Laboratory, UCLA
% Copyright UCLA SPAPL 2009

% first, get all maxima
ymin = min(y);
y = y - ymin;
dy = [y(1); transpose(diff(y))];   % diff() gives a shift by 1: insert one sample
inx1 = diff(dy < 0);    % is +1 for maxima, -1 for minima
inx1(1) = 0; % do not accept maxima at begining
inx1(end) = 0; % do not accept maxima at end
inx = inx1 > 0;         % choose maxima only
ymax= y(inx);
inx = find(inx);
%plot(y);
%hold on;
%plot(inx,ymax,'or');
%hold off;

nofmax = length(ymax);
if nofmax==1
   return;
end
% now filter maxima with window of length winlen
for cnt = 1 : nofmax
   arr = inx(1:cnt);
   cmp = inx(cnt)-winlen;
   arr2 = arr>cmp;
   %ymax(arr2)
   [m, mi] = max(ymax(arr2));
   ymax(arr2)=-60000;
   ymax(mi+length(arr2)-sum(arr2))=m;
   %ymax(arr2)
end
temp = find(ymax>0);
inx  = inx(temp);
ymax = ymax(temp);
%plot(y);
%hold on;
%plot(inx,ymax,'or');
%hold off;
ymax = ymax + ymin;

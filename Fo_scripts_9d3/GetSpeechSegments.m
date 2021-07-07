function [isegmentstart,ksegmentstart]=GetSpeechSegments(...
   kspeechstart,ispeechstart,sig_dB,t2,seg_break)
% GetSpeechSegments- finds speech segments for a phrase or breath based the
% size of the pauses.
%
% Needs the speech area markers kspeechstart, ispeechstart, t2 (time
% corresponding to the kspeechstart), seg_break (size of pause to look for
% in msec)
%
% Returns isegmentstart (index of segments start and end), ksegmentstart
% (a logical array the length of t2 where 1's indicate a start or end of a
% speech segment)
%
% Eric Hunter, NCVS@DCPA, 20060418

RMS_window=mean(diff(t2));
% mark voicing segments 
ksegmentstart=kspeechstart;
if length(ispeechstart)>3  %remove small silences
   for n=2:2:length(ispeechstart)-1
      %check to see if silence is shorter than segment break threshold
%       if ((ispeechstart(n+1)-ispeechstart(n))/Fs)<seg_break/1000
      if ((ispeechstart(n+1)-ispeechstart(n))*RMS_window)<seg_break/1000         
         ksegmentstart(ispeechstart(n))=0;   %remove the silence start
         ksegmentstart(ispeechstart(n+1))=0;   %remove the silence end
      end
   end
end
isegmentstart = find(ksegmentstart == 1);
if length(isegmentstart)>3  %remove small speech
   for n=1:2:length(isegmentstart)
      %check to see if speech is shorter than segment break threshold
      if ((isegmentstart(n+1)-isegmentstart(n))*RMS_window)<seg_break/1000/6         
         ksegmentstart(isegmentstart(n+1))=0;   %remove the silence start
         ksegmentstart(isegmentstart(n))=0;   %remove the silence end
      end
   end
end
isegmentstart = find(ksegmentstart == 1);

if sum(ksegmentstart)==0    % if non, then mark whole segment
   isegmentstart=[1 length(kspeechstart)];
   ksegmentstart=0.*kspeechstart;
   ksegmentstart(1)=1;ksegmentstart(end)=1;
end
if length(isegmentstart)<2     %if only one, make two but mark section w/speech
   if sum(kspeechstart(1:isegmentstart(1)))>...
         sum(kspeechstart(isegmentstart(1):end))  %which has voicing
      isegmentstart=[1 isegmentstart(1)];
      ksegmentstart=[1 ksegmentstart(1)];
   else
      isegmentstart=[isegmentstart(1) length(isegmentstart)];
      ksegmentstart=[ksegmentstart(1) length(ksegmentstart)];
   end  
end
if length(isegmentstart)==2     % if only two, make sure voice is inbetween
   if sum(kspeechstart(1:isegmentstart(1)))>...
         sum(kspeechstart(isegmentstart(1):isegmentstart(2)-1))  %which has voicing
      disp('ODD SPEECH FILE!!!!  YOU SHOULD CHECK BEFORE TRUSTING ANALYSIS')
%       error(' BAD SPEECH FILE!!!!  can''t analyze well, check file in sound editor')
      isegmentstart=[1 isegmentstart(1) isegmentstart(2) length(isegmentstart)];
      ksegmentstart=[1 ksegmentstart(1) ksegmentstart(2) length(ksegmentstart)];
   end
end

if length(isegmentstart)>2     % if at least one segment, check length and such
   if mod(sum(ksegmentstart),2)==1, % if odd, then check to see if the first section has voicing
      if mean(sig_dB(isegmentstart(1):isegmentstart(2)-1))>...
            mean(sig_dB(isegmentstart(2):isegmentstart(3)-1))  %which has voicing
         ksegmentstart(isegmentstart(end))=0;
         isegmentstart=isegmentstart(1:end-1);
      else
         ksegmentstart(isegmentstart(1))=0;
         isegmentstart=isegmentstart(2:end);
      end
   else
      if mean(sig_dB(isegmentstart(1):isegmentstart(2)-1))<...
            mean(sig_dB(isegmentstart(2):isegmentstart(3)-1))  %which has voicing
         ksegmentstart(isegmentstart(end))=0;
         ksegmentstart(isegmentstart(1))=0;
         isegmentstart=isegmentstart(2:end-1);
      end
   end
end
function [tmpstatsHz,tmpstatsST] = freq_stats(fo_array,time_array,freq_Ref)
%freq_stats -- brings in a frequency array and its associated time array and
%calculates statistics both in hz and in semitones
%
%  fo_array - array of values in hz
%  time_array - associate fo array with time
%  freq_Ref - semitone reference

fo_array=fo_array(~isnan(fo_array)); % remove all NaN values
fo_array=fo_array(fo_array>0); % remove all NaN values

tmpstatsHz= basicstats(fo_array);
fo_arrayST = 12*log2(double(fo_array)./freq_Ref); 
tmpstatsST=basicstats(fo_arrayST); % semitones

if ~isnan(time_array)
   time_array=time_array(~isnan(fo_array)); %removes all NaN values
   p = polyfit(time_array,fo_array,1);
   tmpstatsHz.slope=p;
   p = polyfit(time_array,fo_arrayST,1);
   tmpstatsST.slope=p;
else
   tmpstatsHz.slope=NaN;
   tmpstatsST.slope=NaN;
end

tmp=freq_Ref*2.^([tmpstatsST.mean tmpstatsST.mean-tmpstatsST.std tmpstatsST.mean+tmpstatsST.std...
   tmpstatsST.mean-2*tmpstatsST.std tmpstatsST.mean+2*tmpstatsST.std ...
   tmpstatsST.mean-3*tmpstatsST.std tmpstatsST.mean+3*tmpstatsST.std]/12);
tmpstatsHz.stmean=tmp(1);  %mean hz from st calculated stats
tmpstatsHz.stmeanmstd=tmp(2);  %mean-std in semitones then to hz
tmpstatsHz.stmeanpstd=tmp(3);  %mean+std in semitones then to hz
tmpstatsHz.stmeanm2std=tmp(4);  %mean-std in semitones then to hz
tmpstatsHz.stmeanp2std=tmp(5);  %mean+std in semitones then to hz
tmpstatsHz.stmeanm3std=tmp(6);  %mean-std in semitones then to hz
tmpstatsHz.stmeanp3std=tmp(7);  %mean+std in semitones then to hz

end


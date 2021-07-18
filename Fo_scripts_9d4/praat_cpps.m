function [outavqi] = praat_cpps(siga, Fs, fname, code, files)
% function gets the cpps from pratt

%% prepare filename
if isunix == 1; tmp=strrep(cd,'/','.');
elseif ispc == 1; tmp=strrep(cd,'\','.');
end
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav=tmp(3:end);
if length(tmppraatwav)>25
    tmppraatwav=tmppraatwav(length(tmppraatwav)-25:end);
end
uniqID=num2str(randi([10000000 40000000],1,1));
tmppraatwav=['1avqi---' uniqID '---' tmppraatwav '---' '.wav'];


%% create file
audiowrite(tmppraatwav,siga/max(abs(siga)),Fs)
if isunix == 1; movefile([files '/' tmppraatwav],[code '/' tmppraatwav]);
elseif ispc == 1; movefile([files '\' tmppraatwav],[code '\' tmppraatwav]);
end
%% cpps

[s, w] =  system(['cd ' code ' & Praat.exe --run praat_avqi2.praat ' tmppraatwav]);
if isunix == 1; delete([code '/' tmppraatwav])
elseif ispc == 1; delete([code '\' tmppraatwav])
end


%%
if s == -1
    disp('------ AVQI unsuccessful ------')
    outavqi.cpps = NaN;
    outavqi.hnr = NaN;
    outavqi.shim = NaN;
    outavqi.shdb = NaN;
    outavqi.slope = NaN;
    outavqi.tilt = NaN;
    outavqi.avqi = NaN;
else
    datatemp = importdata([code '\avqi2.txt']);
    outavqi.cpps = datatemp(1);
    outavqi.hnr = datatemp(2);
    outavqi.shim = datatemp(3);
    outavqi.shdb = datatemp(4);
    outavqi.slope = datatemp(5);
    outavqi.tilt = datatemp(6);
    outavqi.avqi = NaN;
    delete([code '\avqi2.txt'])
end

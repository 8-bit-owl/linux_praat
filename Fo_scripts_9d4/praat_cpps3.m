function [outavqi] = praat_cpps3(siga,varargin)%Fs, praat_script, code)
% function gets the cpps from pratt

%%
cd0=cd;     %keep track of original folder

%% prepare filename
if isunix == 1; tmp=strrep(cd,'/','.');
elseif ispc == 1; tmp=strrep(cd,'\','.');
end
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav=tmp(3:end);
if length(tmppraatwav)>12
    tmppraatwav=tmppraatwav(length(tmppraatwav)-12:end);
end
uniqID=num2str(randi([10000000 40000000],1,1));
tmppraatwav=['1avqi-' uniqID '-' tmppraatwav '-' ];

Fs=varargin{1};
praat_script=varargin{2};
code= (varargin{3});
foS_lower=varargin{4};
foS_upper=varargin{5};

if ispc == 1; cd( code)
end
%cd('./')

%% prepare for analysis, write out temporary wave file

if (~isa(siga, 'double'))
    error('First input argument must be a vector of doubles');
end
siga(siga == 1) = 32767/32768; %make anythint that is 1 to just less than 1
siga(siga == -1) = -32767/32768; %make anythint that is 1 to just less than 1

audiowrite([tmppraatwav '.wav'], siga, Fs);
% movefile([files '\' tmppraatwav],[code '\' tmppraatwav]);

%% cpps
commstring = ['Praat.exe --run ' praat_script ' '  tmppraatwav '.wav ' tmppraatwav '.txt '...
     ' ' num2str(foS_lower) ' ' num2str(foS_upper) ];

[s, w] =  system(['cd ' code ' & ' commstring]);

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
    if isunix == 1; datatemp = importdata([ code '/' tmppraatwav '.txt ']);
	elseif ispc == 1; datatemp = importdata([code '\' tmppraatwav '.txt ']);
	end
    outavqi.cpps = datatemp(1);
    outavqi.hnr = datatemp(2);
    outavqi.shim = datatemp(3);
    outavqi.shdb = datatemp(4);
    outavqi.slope = datatemp(5);
    outavqi.tilt = datatemp(6);
    outavqi.avqi = NaN;
%     delete([code '\avqi2.txt'])
end

%% clean up the temporary files
delete([tmppraatwav '.txt']);
delete([tmppraatwav '.wav']);
cd(cd0)

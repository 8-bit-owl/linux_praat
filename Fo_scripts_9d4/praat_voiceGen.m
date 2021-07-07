function [perturb] = praat_voiceGen(data, varargin)
% Invokes Praat's perturbative voice analysis methods for a given speech
% signal. Uses the default algorithm parameter values provided with Praat.
%
% Usage:
% [t, F0] = praat_pitch(x, fs)
% [t, F0] = praat_pitch(fn)
% Inputs
%    data      - input signal: must be a row vector
%    varargin(1) - input signal sample rate
%    varargin(2) - where to put the file and do the anlaysis
%    varargin(3) - praat script to write
%    varargin(4) - fo_lower
%    varargin(5) - fo_upper
%    varargin(6) - analysistype  (usually ac, but cc is useful for fry)
%
% Outputs:
%    perturb - structure containing perturbation voice measures
%
% updated by Eric Hunter multiple times

%%
cd0=cd;     %keep track of original folder

%% prepare parameters
Fs=varargin{1};
cd( varargin{2})
praat_script=varargin{3};
fo_lower=varargin{4};
fo_upper=varargin{5};
analysistype=varargin{6};
VoiceThresh=varargin{7};
% get make a file name
tmp=strrep(cd,'/','.');
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav=tmp(3:end); clear tmp
if length(tmppraatwav)>25 % part of filename should be the directory
    tmppraatwav=tmppraatwav(length(tmppraatwav)-25:end);
end

% test if file exist
tmp=dir;
found1=1;
while found1==1
    uniqID=num2str(randi([50000000 90000000],1,1));
    tmppraatwav=['1voice---' uniqID '---' tmppraatwav '---' ];
    found1=0;
for n=1:length(tmp)
    if strcmp(tmp(n).name,[tmppraatwav '_' praat_script])==1
       found1=1;
       uniqID=num2str(randi([50000000 90000000],1,1));
       tmppraatwav=['1voice---' uniqID '---' tmppraatwav '---' ];
       n=length(tmp);
    end
end
end
clear tmp 

%% make a unique praat file, load in generic and replace components

fid = fopen(praat_script,'r');
n = 1; clear A
tline = fgetl(fid);
A{n} = tline;
while ischar(tline)
    n = n+1;
    tline = fgetl(fid);
    A{n} = tline;
end
fclose(fid);

% disp(' ')
% for n=1:numel(A)-1
%     disp(A{n})
% end
% Change cell A
A{7}=strrep(A{7},'fo_lower',num2str( fo_lower));
A{7}=strrep(A{7},'fo_upper',num2str( fo_upper));
A{7}=strrep(A{7},'VoiceThresh',num2str( VoiceThresh));
A{7}=strrep(A{7},'analysistype',( analysistype));
A{9}=strrep(A{9},'fo_lower',num2str( fo_lower));
A{9}=strrep(A{9},'fo_upper',num2str( fo_upper));
A{9}=strrep(A{9},'analysistype',( analysistype));
A{14}=strrep(A{14},'fo_lower',num2str( fo_lower));
A{14}=strrep(A{14},'fo_upper',num2str( fo_upper));
A{14}=strrep(A{14},'VoiceThresh',num2str( VoiceThresh));
A{14}=strrep(A{14},'analysistype',( analysistype));
for n=32:45
A{n}=strrep(A{n},'voice.txt',[tmppraatwav '.txt']);
end
% disp(' ')
% for n=1:numel(A)-1
%     disp(A{n})
% end


%% write out new file
tmppraatwav_psc=[tmppraatwav '_' praat_script];
fid = fopen(tmppraatwav_psc, 'w');
for n = 1:numel(A)
    if A{n+1} == -1
        fprintf(fid,'%s', A{n});
        break
    else
        fprintf(fid,'%s\n', A{n});
    end
end
fclose(fid);
% open(tmppraatwav_psc)

%% prepare for analysis, write out temporary wave file

if (~isa(data, 'double'))
    error('First input argument must be a vector of doubles');
end
data(data == 1) = 32767/32768; %make anythint that is 1 to just less than 1
data(data == -1) = -32767/32768; %make anythint that is 1 to just less than 1

audiowrite([tmppraatwav '.wav'], data, Fs);

%% analysis
commstring = ['praatcon ' tmppraatwav_psc ' ' tmppraatwav '.wav'];

%     commstring = 'praatcon_win98 praat_pitch.psc praat.wav';
[s, w] = system(commstring,'-echo');
 
%% load the results
% results = load([tmppraatwav '.txt']);

fieldnames = {'jitter','jitter_abs','jitter_rap','jitter_ppq5','jitter_ddp',...
    'shimmer','shimmer_db','shimmer_apq3','shimmer_apq5','shimmer_apq11','shimmer_dda',...
    'nhr','hnr'};

fh = fopen([tmppraatwav '.txt']);
perturb = [];
for j = 1:length(fieldnames)
    line = fgetl(fh);
    if (strmatch(line, '--undefined--'))
        perturb = setfield(perturb, fieldnames{j}, NaN);
    else
        perturb = setfield(perturb, fieldnames{j}, str2num(line));
    end
end
fclose(fh);

%% clean up the temporary files
delete(tmppraatwav_psc);
delete([tmppraatwav '.txt']);
delete([tmppraatwav '.wav']);

cd(cd0)

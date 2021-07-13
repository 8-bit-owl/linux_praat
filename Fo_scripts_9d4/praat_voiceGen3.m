function [perturb] = praat_voiceGen3(data, varargin)
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

%% get & make a file name
tmp=strrep(cd,'.','/')
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav0=tmp(3:end); clear tmp
if length(tmppraatwav0)>12 % part of filename should be the directory
    tmppraatwav0=tmppraatwav0(length(tmppraatwav0)-12:end);
end

%% prepare parameters
Fs=varargin{1};
cd( varargin{2})
praat_script=varargin{3};
foS_Analysis=varargin{4};
foS_step=varargin{5};
foS_lower=varargin{6};
foS_upper=varargin{7};
foS_NCand=varargin{8};
foS_Accuracy=varargin{9};
foS_SilenceThrsh=varargin{10};
foS_VoiceThrsh=varargin{11};
foS_OctCost=varargin{12};
foS_OctJumpCost=varargin{13};
foS_VoiceUnvoiceCost=varargin{14};
pertS_MaxPeriodFact=varargin{15};
pertS_MaxAmpFact=varargin{16};

% voiceReport$ = Voice report... 0 0 fo_lower fo_upper 1.3 1.6 0.03 fo_VoiceThrsh

%% test if file exist
tmp=dir;
found1=1;
while found1==1
    uniqID=num2str(randi([5000000000 9000000000],1,1));
    tmppraatwav=['1pert-' uniqID '-' tmppraatwav0 '-' ];
    found1=0;
for n=1:length(tmp)
    if strcmp(tmp(n).name,[tmppraatwav '_' praat_script])==1
       found1=1;
       uniqID=num2str(randi([5000000000 9000000000],1,1));
       tmppraatwav=['1pert-' uniqID '-' tmppraatwav0 '-' ];
       n=length(tmp);
    end
end
end
clear tmp 


%% prepare for analysis, write out temporary wave file
if (~isa(data, 'double'))
    error('First input argument must be a vector of doubles');
end
data(data == 1) = 32767/32768; %make anythint that is 1 to just less than 1
data(data == -1) = -32767/32768; %make anythint that is 1 to just less than 1

audiowrite([tmppraatwav '.wav'], data, Fs);

%% analysis
% commstring = ['praatcon ' tmppraatwav_psc ' ' tmppraatwav '.wav'];

commstring = ['praatcon ' praat_script ' '  tmppraatwav '.wav ' tmppraatwav '.txt ' ...
    num2str(foS_lower) ' ' num2str(foS_upper) ' ' num2str(foS_NCand) ' ' ...
    foS_Accuracy ' ' num2str(foS_SilenceThrsh) ' '  num2str(foS_VoiceThrsh) ' ' ...
    num2str(foS_OctCost) ' ' num2str(foS_OctJumpCost) ' '  num2str(foS_VoiceUnvoiceCost) ' ' ...
    num2str(foS_step) ' ' foS_Analysis ' ' num2str(pertS_MaxPeriodFact) ' ' num2str(pertS_MaxAmpFact)];

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
delete([tmppraatwav '.txt']);
delete([tmppraatwav '.wav']);

cd(cd0)

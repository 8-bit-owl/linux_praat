function [t, F0] = praat_pitchGen3(data, varargin)
% Invokes Praat's pitch (F0) extraction method for a given speech signal.
% Uses the default pitch extraction algorithm parameter values provided
% with Praat.
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
%    varargin(7) - VoiceThresh usually 0.45
%    http://www.fon.hum.uva.nl/praat/manual/Sound__To_Pitch__ac____.html
%
% Outputs:
%    t      - vector of pitch mark time instants
%    F0     - vector of F0 estimates at the pitch mark time instants
%
% updated by Eric Hunter multiple times

%%
cd0=cd;     %keep track of original folder

%% get & make a file name
tmp=strrep(cd,'/','.');
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav0=tmp(3:end);
if length(tmppraatwav0)>12
    tmppraatwav0=tmppraatwav0(length(tmppraatwav0)-12:end);
end

%% prepare Parameters
Fs=varargin{1};
persistent x_variable;
if isempty(x_variable)
    cd ( varargin{2})
    x_variable = 0;
end
%cd ( varargin{2}) 
praat_script=varargin{3};
pitchAnalysis=varargin{4};
fo_step=varargin{5};
fo_lower=varargin{6};
fo_upper=varargin{7};
pitchNCand=varargin{8};
pitchAccuracy=varargin{9};
pitchSilenceThrsh=varargin{10};
pitchVoiceThrsh=varargin{11};
pitchOctCost=varargin{12};
pitchOctJumpCost=varargin{13};
pitchVoiceUnvoiceCost=varargin{14};

%% test if file exist
tmp=dir;
found1=1;
while found1==1
    uniqID=num2str(randi([10000000000 40000000000],1,1));
    tmppraatwav=['1fo-' uniqID '-dir-' tmppraatwav0 '-' ];
    found1=0;
for n=1:length(tmp)
    if strcmp(tmp(n).name,[tmppraatwav '_' praat_script])==1
       found1=1;
       uniqID=num2str(randi([10000000000 40000000000],1,1));
       tmppraatwav=['1voice---' uniqID '---' tmppraatwav0 '---' ];
       n=length(tmp);
    end
end
end
clear tmp uniqID


%% prepare for analysis, write out temporary wave file

if (~isa(data, 'double'))
    error('First input argument must be a vector of doubles');
end
data(data == 1) = 32767/32768; %make anythint that is 1 to just less than 1
data(data == -1) = -32767/32768; %make anythint that is 1 to just less than 1

audiowrite([tmppraatwav '.wav'], data, Fs);

%% analysis

commstring = ['./praat_barren ' praat_script ' '  tmppraatwav '.wav ' tmppraatwav '.txt ' ...
    num2str(fo_lower) ' ' num2str(fo_upper) ' ' num2str(pitchNCand) ' ' ...
    num2str(pitchAccuracy) ' ' num2str(pitchSilenceThrsh) ' '  num2str(pitchVoiceThrsh) ' ' ...
    num2str(pitchOctCost) ' ' num2str(pitchOctJumpCost) ' '  num2str(pitchVoiceUnvoiceCost) ' ' ...
    num2str(fo_step) ' ' pitchAnalysis];

%     commstring = 'praatcon_win98 praat_pitch.psc praat.wav';
[s, w] = system(commstring,'-echo');
 
results = load([tmppraatwav '.txt']);

%% clean up the temporary files
delete([tmppraatwav '.txt']);
delete([tmppraatwav '.wav']);

%% adjust the results and return to original directory
if (~isempty(results))
    t  = results(:,1);
    F0 = results(:,2);
else
    warning('No pitch mark data obtained from Praat');
    t  = NaN;
    F0 = NaN;
end

cd(cd0)

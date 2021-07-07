function [t, F0] = praat_pitchGen(data, varargin)
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

%% prepare parameters
tmp=strrep(cd,'/','.');
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav=tmp(3:end);
if length(tmppraatwav)>25
    tmppraatwav=tmppraatwav(length(tmppraatwav)-25:end);
end
uniqID=num2str(randi([10000000 40000000],1,1));
tmppraatwav=['1pitch---' uniqID '---' tmppraatwav '---' ];

Fs=varargin{1};
cd( varargin{2})
praat_script=varargin{3};
fo_lower=varargin{4};
fo_upper=varargin{5};
analysistype=varargin{6};
VoiceThresh=varargin{7};

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
A{5}=strrep(A{5},'fo_lower',num2str( fo_lower));
A{5}=strrep(A{5},'fo_upper',num2str( fo_upper));
A{5}=strrep(A{5},'analysistype',num2str( analysistype));
A{5}=strrep(A{5},'VoiceThresh',num2str( VoiceThresh));
A{7}=strrep(A{7},'pitch.txt',[tmppraatwav '.txt']);
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
 
results = load([tmppraatwav '.txt']);

%% clean up the temporary files
delete(tmppraatwav_psc);
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

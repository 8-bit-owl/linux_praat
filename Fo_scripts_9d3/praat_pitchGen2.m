function [t, F0] = praat_pitchGen2(data, varargin)
% Invokes Praat's pitch (F0) extraction method for a given speech signal.
% Uses the default pitch extraction algorithm parameter values provided
% with Praat.
%
% Usage:
% [t, F0] = praat_pitch(x, fs)
% [t, F0] = praat_pitch(fn)
% Inputs
%    data      - input signal: must be a row vector
%    Fs=varargin{1};
%    cd(varargin{2})
%    praat_script=varargin{3};
%    pitchAnalysis=varargin{4};
%    fo_step=varargin{5};
%    fo_lower=varargin{6};
%    fo_upper=varargin{7};
%    pitchNCand=varargin{8};
%    pitchAccuracy=varargin{9};
%    pitchSilenceThrsh=varargin{10};
%    pitchVoiceThrsh=varargin{11};
%    pitchOctCost=varargin{12};
%    pitchOctJumpCost=varargin{13};
%    pitchVoiceUnvoiceCost=varargin{14};
%
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
tmp=strrep(cd,'\','.');
tmp=strrep(tmp,' ','_');
tmp=strrep(tmp,'_','-');
tmppraatwav=tmp(3:end);
if length(tmppraatwav)>25
    tmppraatwav=tmppraatwav(length(tmppraatwav)-25:end);
end
uniqID=num2str(randi([10000000 40000000],1,1));
tmppraatwav=['1pitch---' uniqID '---' tmppraatwav '---' ];

Fs=varargin{1};
cd(varargin{2})
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
A{5}=strrep(A{5},'fo_step',num2str( fo_step));
A{5}=strrep(A{5},'fo_lower',num2str( fo_lower));
A{5}=strrep(A{5},'fo_upper',num2str( fo_upper));
A{5}=strrep(A{5},'pitchAnalysis',pitchAnalysis);
A{5}=strrep(A{5},'pitchVoiceThrsh',num2str( pitchVoiceThrsh));
A{5}=strrep(A{5},'pitchNCand',num2str( pitchNCand));
A{5}=strrep(A{5},'pitchAccuracy',num2str( pitchAccuracy));
A{5}=strrep(A{5},'pitchSilenceThrsh',num2str( pitchSilenceThrsh));
A{5}=strrep(A{5},'pitchVoiceThrsh',num2str( pitchVoiceThrsh));
A{5}=strrep(A{5},'pitchOctCost',num2str( pitchOctCost));
A{5}=strrep(A{5},'pitchOctJumpCost',num2str( pitchOctJumpCost));
A{5}=strrep(A{5},'pitchVoiceUnvoiceCost',num2str( pitchVoiceUnvoiceCost));
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
[s, w] = dos(commstring,'-echo');
 
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

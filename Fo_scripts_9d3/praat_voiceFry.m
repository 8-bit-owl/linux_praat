function perturb = praat_voiceFry(data, varargin)
% Invokes Praat's perturbative voice analysis methods for a given speech
% signal. Uses the default algorithm parameter values provided with Praat.
%
% Usage:
% perturb = praat_voice(x, fs)
% perturb = praat_voice(fn)
% Inputs
%    x       - input signal: must be a row vector
%    fs      - input signal sample rate
%    fn      - WAV filename
%
% Outputs:
%    perturb - structure containing perturbation voice measures
%
% (c) 2008 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, E.J. Hunter, J. Spielman, L.O. Ramig (2008)
% Suitability of dysphonia measurements for telemonitoring of Parkinson’s
% disease
% IEEE Transactions on Biomedical Engineering
%
% added other metrics.  20130607  - Eric Hunter

if (nargin == 2)
    if (~isa(data, 'double'))
        error('First input argument must be a vector of doubles');
    end
    data(data == 1) = 32767/32768;
%     wavwrite(data, varargin{1}, 16, 'praat.wav');
    audiowrite('praat.wav',data, varargin{1});
    commstring = 'praatcon praat_voiceFry.psc praat.wav';
%     commstring = 'praatcon_win98 praat_voiceF.psc praat.wav';    
    
[s, w] = dos(commstring);
delete('praat.wav');
elseif (nargin == 3)
    if (~isa(data, 'double'))
        error('First input argument must be a vector of doubles');
    end
    data(data == 1) = 32767/32768;
    %     wavwrite(data, varargin{1}, 16, 'praat.wav');
    audiowrite([varargin{2} '\praat.wav'],data, varargin{1});
    commstring = 'praatcon praat_voiceFry.psc praat.wav';
    %     commstring = 'praatcon_win98 praat_pitch.psc praat.wav';
    [s, w] = dos(['cd ' varargin{2} ' & ' commstring]);
    delete([varargin{2} '\praat.wav']);
elseif (nargin == 1)
    if (isa(data,'char'))
        commstring = ['praatcon praat_voiceFry.psc ' data];
%         commstring = ['praatcon_win98 praat_voiceF.psc ' data];        
    else
        error('First input argument must be a filename');
    end
    [s, w] = dos(commstring);
else
    error('Incorrect number of arguments');
end

fieldnames = {'jitter','jitter_abs','jitter_rap','jitter_ppq5','jitter_ddp',...
              'shimmer','shimmer_db','shimmer_apq3','shimmer_apq5','shimmer_apq11','shimmer_dda',...
              'nhr','hnr'};

if (nargin < 3)
    fh = fopen('voice.txt');
else
    fh = fopen([varargin{2} '\voice.txt']);
end
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
if (nargin < 3)
    delete('voice.txt');
else
    delete([varargin{2} '\voice.txt']);
end

function [t, F0] = praat_pitchFry(data, varargin)
% Invokes Praat's pitch (F0) extraction method for a given speech signal.
% Uses the default pitch extraction algorithm parameter values provided
% with Praat.
%
% Usage:
% [t, F0] = praat_pitch(x, fs)
% [t, F0] = praat_pitch(fn)
% Inputs
%    x      - input signal: must be a row vector
%    fs     - input signal sample rate
%    fn     - WAV filename
%
% Outputs:
%    t      - vector of pitch mark time instants
%    F0     - vector of F0 estimates at the pitch mark time instants
%
% (c) 2008 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, E.J. Hunter, J. Spielman, L.O. Ramig (2008)
% Suitability of dysphonia measurements for telemonitoring of Parkinson’s
% disease
% IEEE Transactions on Biomedical Engineering

%%

if (nargin == 2)
    if (~isa(data, 'double'))
        error('First input argument must be a vector of doubles');
    end
    data(data == 1) = 32767/32768;
%     wavwrite(data, varargin{1}, 16, 'praat.wav');
    audiowrite('praat.wav',data, varargin{1});
    commstring = 'praatcon praat_pitchFry.psc praat.wav';
%     commstring = 'praatcon_win98 praat_pitch.psc praat.wav';    
    [s, w] = dos(commstring);
    delete('praat.wav');
elseif (nargin == 3)
    if (~isa(data, 'double'))
        error('First input argument must be a vector of doubles');
    end
    data(data == 1) = 32767/32768;
%     wavwrite(data, varargin{1}, 16, 'praat.wav');
    audiowrite([varargin{2} '\praat.wav'],data, varargin{1});
    commstring = 'praat praat_pitchFry.psc praat.wav';
%     commstring = 'praatcon_win98 praat_pitch.psc praat.wav';    
    [s, w] = dos(['cd ' varargin{2} ' & ' commstring],'-echo');
    delete([varargin{2} '\praat.wav']);
elseif (nargin == 1)
    if (isa(data,'char'))
        commstring = ['praatcon praat_pitchFry.psc ' data];
        commstring = ['praatcon_win98 praat_pitchFry.psc ' data];
    else
        error('First input argument must be a filename');
    end
    [s, w] = dos(commstring);
else
    error('Incorrect number of arguments');
end

if nargin < 3
    results = load('pitch.txt');
    delete('pitch.txt');
else
    results = load([varargin{2} '\pitch.txt']);
    delete([varargin{2} '\pitch.txt']);
end

if (~isempty(results))
    t  = results(:,1);
    F0 = results(:,2);
else
    warning('No pitch mark data obtained from Praat');
    t  = NaN;
    F0 = NaN;
end

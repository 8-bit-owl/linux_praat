function measures = voice_measures3(x, fs, varargin)
% Calculates a range of classical and non-classical voice measures for a
% given input voice recording, using standard algorithm parameter values.
% several of the commands come from here:
% http://www.fon.hum.uva.nl/praat/manual/Voice_2__Jitter.html
% Usage:
% measures = voice_measures(x, fs)
% Inputs
%    x        - input signal: must be a row vector
%    fs       - input signal sample rate
%
% Outputs:
%    measures - structure containing obtained voice measures
%
% (c) 2008 Max Little. If you use this code, please cite:
% M.A. Little, P.E. McSharry, E.J. Hunter, J. Spielman, L.O. Ramig (2008)
% Suitability of dysphonia measurements for telemonitoring of Parkinson’s
% disease
% IEEE Transactions on Biomedical Engineering
%
% added other metrics.  20130607  - Eric Hunter
% Full Rewrite.  20180907  - Eric Hunter

% Default voice analysis parameter values at predefined sample rate
Tmax     = 1000;
Fs       = 25000;   % Sample rate
d        = 4;       % Embedding dimension
tau      = 50;      % Embedding delay
eta      = 0.2;     % RPDE close returns radius
dfa_scaling = (50:20:200)';   % DFA scaling range
Ssample  = 5000;    % Start sample
Esample  = 20000;   % End sample

% Have to guard against silent and zero-sized files
if (length(x) > Esample)
 
    % Resample to analysis rate
    x = resample(x, Fs, fs);
    x = x(Ssample:min([Esample length(x)]));

    % Amplitude normalization
    x = x - mean(x);
    x = x/max(abs(x));
   
    
    % Get RPDE and DFA
    try
        H_norm = rpde(x, d, tau, eta, Tmax);
        [dfa, iv, fl] = fastdfa(x, dfa_scaling);
        alpha_norm    = 1/(1+exp(-dfa));
    catch
        disp('skipping RPD and DFA')
        H_norm = NaN;
        dfa=NaN;
        alpha_norm=NaN;
    end
    
    % Get PPE
    try
        E_norm = ppe(x, Fs);
    catch
        E_norm=NaN;
        disp('skipping PPE')
    end
    
    % Prepare output
    measures = perturb;
    measures.AlphaRatio = AlphaRatio;
%     measures.LTASf = LTAS.f;
%     measures.LTASdBspec = LTAS.dBspectrum;
    measures.LTASdBfft=LTAS.totdBfft;
    measures.LTASnFFT=LTAS.Nfft;
    measures.LTASframes=LTAS.frames;
    measures.LTASdBspect_slope=LTAS.dBspect_slope;
    
    measures.rpde = H_norm;
    measures.dfa = alpha_norm;
    measures.ppe = E_norm;
else
    warning('Signal length not long enough to calculate measures');
    measures.jitter = NaN;
end

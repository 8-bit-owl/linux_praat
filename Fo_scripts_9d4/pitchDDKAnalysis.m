function [results,parameters] = pitchDDKAnalysis(x,fs,parameters)
% This function analyzes pitch DDK recordings.
% Input:
%    x -- Nx1 real vector, time domain signal
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of parameters
% Output:
%    results -- struct of results

% Set default parameters:
if nargin<3,
   parameters(1).plim = [75 500]; % [min max] Hz
   parameters(1).dt = 0.01; % sec
   parameters(1).dlog2p = 1/48; % 1/48-octave stepsize
   parameters(1).dERBs = 1/20; % ERBrate
   parameters(1).woverlap = 0.5; % window overlap fraction
   parameters(1).sTHR = 0.2; % pitch strength threshold
   parameters(1).plot = true; % plot
end;

% Resample to 20 kHz:
if fs~=20e3,
   % Get integer ratio of 20e3/fs:
   [a,b] = rat(20e3/fs); % a/b ~ 20e3/fs
   
   % Resample:
   x = resample(x,a,b);
   fs = 20e3; % Hz;
end;

% Get pitch track:
[P,T] = audswipep(x,fs,parameters.plim,parameters.dt,parameters.dlog2p,parameters.dERBs,...
   parameters.woverlap,parameters.sTHR); % P Hz, T sec

% Convert P from Hz to semitones:
C0 = 440*2^(-57/12); % C0 is 57 semitones below A4 (440 Hz), C0~16.3516 Hz
P = 12*log(P/C0)/log(2); % C0 reference

% Find slope of P:
dP = diff(P)/(T(2)-T(1)); % semitone/sec

% Calculate moments of dP:
dP1 = dP(~isnan(dP));
dPmean = mean(dP1);
dPvar = var(dP1);
dPstd = std(dP1);
dPskew = skewness(dP1);
dPkurt = kurtosis(dP1);

% Calculate maxfilter of dP:
dPmaxfilt = maxfilter(dP,80);

% Save to output:
results(1).dPmean = dPmean;
results(1).dPvar = dPvar;
results(1).dPstd = dPstd;
results(1).dPskew = dPskew;
results(1).dPkurt = dPkurt;

% Find histogram of dP:
L = length(dP);
Lhist = max(10,min(5000,L/5)); % capped between 10 and 50000
[HdP,bindP] = hist(dP,Lhist);

% Plot results:
if parameters.plot,
   figure; % new figure
   clf;
   subplot(3,1,1); % pitch track
   plot(T,P,'m.-');
   grid on;
   xlabel('Time, sec');
   ylabel('Pitch, semitones');
   subplot(3,1,2); % delta pitch track
   plot(T(2:end),dP,'rd-');
   xlabel('Time, sec');
   ylabel('Delta pitch, semitone/sec');
   grid on;
   hold on
   plot(T(2:end),dPmaxfilt,'g.-');
   subplot(3,1,3); % dP histogram
   plot(bindP,HdP,'ko-');
   xlabel('Delta pitch, semitone/sec');
   ylabel('Histogram');
   grid on;
   title(['Mean: ',num2str(dPmean),', StDev: ',num2str(dPstd),', Skewness: ',num2str(dPskew),...
      ', Kurtosis: ',num2str(dPkurt)]);
end;

return;

function y = maxfilter(x,n)
% This function returns the max of x in windows of length n

% Ensure x is COLUMN vector:
x = x(:);

% Init output:
y = zeros(size(x));

% Get lower/upper window size:
n1 = floor(n/2);
n2 = n-n1-1; % x(n1:n2)

for p=1:length(x),
   y(p) = max(x(max(1,p-n1):min(length(x),p+n2)));
end;

return;

% Bye!
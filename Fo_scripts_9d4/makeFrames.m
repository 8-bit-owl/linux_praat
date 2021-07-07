function [xFrames,t,parameters] = makeFrames(x,fs,parameters)
% This function returns a cell array of analysis frames of x according to parameters. 
% If frameRate = 0, xFrames{1} = x is the only cell in xFrames.
% All frames are of length L = round(fs*windowSize) except the last frame which contains the last
% bit of x with between 1 and L samples.
% The function accounts for fractional frame offset = fs/frameRate (window overlap not a constant
% integer when frame offset is not integer as assumed by functions like spectrogram.m).
% Input [defaults]:
%    x -- Jx1 real vector, time domain signal to analyze, J samples
%    fs -- real scalar, Hz, sampling rate of x
%    parameters -- struct of relevant parameters
%       .windowSize -- real scalar, sec, duration of analysis window [20e-3]
%       .frameRate -- real scalar, frames/sec, number of analysis windows per second [100]
% Output:
%    xFrames -- Tx1 cell array, T analysis frames, L samples of x in each frame except last
%    t -- Tx1 integer vector, index into x of first sample of each frame
%    parameters -- same as input, or defaults if empty

% Set defaults:
defaults = struct([]);
defaults(1).windowSize = 20e-3; % sec
defaults(1).frameRate = 100; % frames/sec

% Check inputs:
if nargin<3,
   parameters = defaults;
end;
x = x(:); % ensure x is COLUMN vector

% Determine number of frames and frame offset:
J = length(x); % samples, data length
L = round(fs*parameters.windowSize); % samples, analysis window length
if parameters.frameRate==0 || L>=J, % one frame only
   T = 1; % frames
   frameOffset = 1; % samples, offset from one frame to the next, arbitrary when T = 1.
else % two or more frames
   % Count number of frames until frame offset pushes last sample of last window past length of x,
   % account for fractional frame offset when fs/frameRate is not integer:
   frameOffset = fs/parameters.frameRate; % samples, offset from one frame to the next, may be fractional
   T = 1; % frames, init
   while L+round(T*frameOffset) < J, % window length + frame offsets of NEXT full frame
      T = T+1; % increment to NEXT full frame
   end;
   T = T+1; % increment to FINAL frame, which may be of length 1 to L
end;

% Separate x into frames:
xFrames = cell(T,1); % init
t = zeros(T,1); % index into x
for p=1:T, % for each analysis frame
   if p<T, % all but final frame
      h = [1:L]+round((p-1)*frameOffset); % samples, indices into x, FULL frame
   else % final frame
      h = [1+round((p-1)*frameOffset):J]; % samples, FULL or PARTIAL frame
   end;
   t(p) = h(1);
   xFrames{p} = x(h); % analysis window
end;

return;

% Bye!
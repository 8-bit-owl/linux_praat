function parameters = promptSWIPEParameters(parameters,defaultsOnly)
% This function prompts for HFCC parameters or returns default values according to input flag.
% Input:
%    parameters -- scalar struct or empty, any input fields
%    defaultsOnly -- logical scalar, [false]: prompt user for parameter values, true: skip prompt
%       and return default values
% Output:
%    parameters -- scalar struct, same as input with added fields

% Check inputs:
if nargin<2,
   defaultsOnly = false; % prompt user for inputs
end;
if nargin<1 || isempty(parameters), % set defaults
   parameters(1).calcSWIPE = true;
   parameters(1).plim = [65 450];
   parameters(1).dt = 1/100; % sec
   parameters(1).dlog2p = 1/24;
   parameters(1).dERBs = 0.1;
   parameters(1).woverlap = 0.5;
   parameters(1).sTHR = 0.1;
end;

% Set display strings for special cases:
calcSWIPEDisplay = 'false';
if parameters.calcSWIPE,
   calcSWIPEDisplay = 'true';
end;
[dlog2pNUM,dlog2pDEN] = rat(parameters.dlog2p);
dlog2pDisplay = [num2str(dlog2pNUM),'/',num2str(dlog2pDEN)]; % display as fraction

% Set parameter prompts and default values:
parameterBlock = {...
   'Calculate SWIPE pitch and pitch strength estimates [true, false]',calcSWIPEDisplay;...
   'MINIMUM pitch frequency range, Hz',num2str(parameters.plim(1));...
   'MAXIMUM pitch frequency range, Hz',num2str(parameters.plim(2));...
   'Frame rate, fps',num2str(1/parameters.dt);...
   'Pitch candidate frequency spacing, octave',dlog2pDisplay;...
   'Spectrum sample frequency spacing, ERBrate',num2str(parameters.dERBs);...
   'Window overlap, range: [0 1]',num2str(parameters.woverlap);...
   'Pitch strength threshold, range: [0 1]',num2str(parameters.sTHR)};
displayString = parameterBlock(:,1);
valueString = parameterBlock(:,2);

% Prompt for parameter settings if necessary:
if ~defaultsOnly, % prompt user for values
   a = inputdlg(displayString,'SWIPE PARAMETERS',[1 80],valueString); % [] if cancel button
else % use default values without prompting
   a = []; % same as cancel button
end;

% Handle cancel button:
if isempty(a), % cancel button
   a = valueString; % use initial values
end;

% Store parameters:
parameters(1).calcSWIPE = str2num(a{1});
parameters(1).plim = [str2num(a{2}) str2num(a{3})];
parameters(1).dt = 1/str2num(a{4}); % sec
parameters(1).dlog2p = str2num(a{5});
parameters(1).dERBs = str2num(a{6});
parameters(1).woverlap = str2num(a{7});
parameters(1).sTHR = str2num(a{8});

return;

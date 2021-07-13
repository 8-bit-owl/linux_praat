function changeFs(rootPath)
% This function prompts for a file or batch of files and prompts for a new sampling rate to apply to
% all files.  The data in the files is resampled up or down to match the new sampling rate.
% Input [default]:
%    rootPath -- string, name of root path to start with in file selection window ['.\']

% Mark Skowronski, February 16, 2014

% Parameters:
defaultFs = 44100; % Hz

% Check input:
if nargin<1,
   rootPath = './'; % current directory, including trailing '\'
end;

% Prompt for file/files to process:
[fileName,pathName] = uigetfile([rootPath,'*.wav'],'Select WAV file(s) for resampling.  Use CTRL or SHIFT keys to select multiple files.','Multiselect','on'); % fileName=string for single-file selection, =cell array for multiple file selection
if ~iscell(fileName), % single file selected or Cancel button
   if fileName(1)==0, % Cancel button
      return; % exit without resampling
   else
      fileName = {fileName}; % Single file selected, convert to cell array for compatibility
   end;
end;

% Prompt for new sampling rate:
a = inputdlg('New sampling rate, Hz','SAMPLING RATE',[1 80],{num2str(defaultFs)}); % [] if cancel button
if isempty(a),
   return; % exit without resampling
else
   fsNew = str2num(a{1});
end;

% Warn that changes are permanent:
buttonLabel = {'Continue','Cancel'};
g = questdlg('WARNING: Resampling permanently changes WAV files.  CONTINUE resampling or CANCEL?',...
   'CONFIRM WAV FILE RESAMPLING',buttonLabel{1},buttonLabel{2},buttonLabel{2});
if isempty(g) || strcmp(g,buttonLabel{2}), % cancel
   return; % exit with resampling
end;

% Resample each selected file, write data to old WAV file, if fs different:
hWait = waitbar(0,'Resampling WAV files'); % init wait bar
L = length(fileName);
for p=1:L, % for each file
   % Update wait bar
   waitbar(p/L,hWait);
   
   % Read WAV file:
   [x,fs] = audioread([pathName,fileName{p}]);
   
   % Resample, save if necessary:
   if fs~=fsNew,
      % Resample:
      [a,b] = rat(fsNew/fs); % smallest integer ratio
      x = resample(x,a,b);
      
      % Write WAV file:
      audiowrite([pathName,fileName{p}],x,fsNew);
   end;
end;

% Close wait bar:
waitbar(1,hWait,['Finished resampling ',num2str(L),' file(s) to ',num2str(fsNew),' Hz.']);
pause(2);
close(hWait); % close wait bar window

return;

% Bye!
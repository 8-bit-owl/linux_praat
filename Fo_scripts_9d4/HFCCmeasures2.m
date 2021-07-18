function output = HFCCmeasures2(parameters)
% This function prompts for a directory of WAV files to process HFCC measures.  The results are
% written to an Excel file in the WAV file directory.
% Input:
%    parameters -- struct of parameters to use for HFCC feature extraction and other sentence
%       processing parameters.  See HFCC.m and processSentence.m for details.

% Check inputs:
if nargin<1,
   parameters = []; % set defaults from processSentence.m and HFCC.m
end;

% Prompt for directory to process:
dirName = uigetdir('.','Select directory of WAV files to process:');
if dirName(1)==0, % user hit Cancel button or closed user interface
   output = []; % init empty output variable
   return;
end;
if isunix == 1; dirName = [dirName,'/']; % append '\' to end of directory name
elseif ispc == 1; dirName = [dirName,'\']; % append '\' to end of directory name
end

% Set output file name:
outputFile = [dirName,'HFCCmeasures.xls'];

% Get list of WAV files:
m = dir([dirName,'*.wav']);

% Process each WAV file:
for p=1:length(m),
   % Read WAV file:
   [x,fs(p)] = wavread([dirName,m(p).name]);
   lenX(p) = length(x);
   
   % Trim endpoints, calculate HFCC measures:
   [output(p),parameters] = processSentences2(m,fs(p),parameters);
end;

% Include file names from m with output:
for p=1:length(m),
   output(p).name = m(p).name;
end;

% Write output to Excel:
header = {'File','Original fs, Hz','Resampled fs, Hz','Original duration, s',...
   'Trimmed duration, s','Active speech level, dB DFS','Activity factor','HFCC st. dev. sum',...
   'Delta HFCC st. dev. sum','Delta Delta HFCC st. dev. sum', 'Filterbank measures', 'Delta FB', 'Delta Delta FB'};
for p=1:length(m),
   t = p+1; % ROW pointer into Excel output cell array
   header{t,1} = upper(m(p).name);
   header{t,2} = fs(p);
   header{t,3} = output(p).fs;
   header{t,4} = lenX(p)/fs(p);
   header{t,5} = output(p).lenXtrimmed;
   header{t,6} = output(p).activeSpeechLevel;
   header{t,7} = output(p).activityFactor;
   header{t,8} = output(p).ccStdSum;
   header{t,9} = output(p).DccStdSum;
   header{t,10} = output(p).DDccStdSum;
   header {t,11} = output(p).FBMeasures;
   header {t,12} = output(p).DFBMeasures;
   header {t,13} = output(p).DDFBMeasures;
end;

% Write output:
xlswrite(outputFile,header,'Sheet1','A1');

% Write HFCC st.dev. for each coefficient index:
header2 = {'File','HFCC st.dev. cc(1)','cc(2)','cc(3)','cc(4)','cc(5)','cc(6)','cc(7)','cc(8)',...
   'cc(9)','cc(10)','cc(11)','cc(12)','cc(13)'};
header3 = {'File','Delta HFCC st.dev. cc(1)','cc(2)','cc(3)','cc(4)','cc(5)','cc(6)','cc(7)','cc(8)',...
   'cc(9)','cc(10)','cc(11)','cc(12)','cc(13)'};
for p=1:length(m),
   t = p+1;
   header2{t,1} = upper(m(p).name);
   header3{t,1} = upper(m(p).name);
   ccStd = std(output(p).HFCC.cc,[],2);
   DccStd = std(output(p).HFCC.Dcc,[],2);
   for p1=1:length(ccStd),
      header2{t,1+p1} = ccStd(p1);
      header3{t,1+p1} = DccStd(p1);
   end;
end;
% Write output:
xlswrite(outputFile,header2,'Sheet2','A1');
xlswrite(outputFile,header3,'Sheet3','A1');

return;

% Bye!
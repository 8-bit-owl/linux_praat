function output = HFCCmeasures3(parameters) 
%Edited by Lisa Kopf on September 14, 2013
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
dirName = [dirName,'\']; % append '\' to end of directory name

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
   [output(p),parameters] = processSentence3(x,fs(p),parameters);
end;

% Include file names from m with output:
for p=1:length(m),
   output(p).name = m(p).name;
end;

% Write output to Excel:
header1 = {'HFCC',' ',' ',' ','Delta HFCC',' ',' ',' ','Delta Delta HFCC'};
header2 = {'File','Original fs, Hz','Resampled fs, Hz','Original duration, s',...
   'Trimmed duration, s','Active speech level, dB DFS','Activity factor','sum SD','sum var.','eucl. dist.','abs. dist.'...
   'sum SD','sum var.','eucl. dist.','abs. dist.','sum SD','sum var.','eucl. dist.','abs. dist.'};
for p=1:length(m),
   t = p+1; % ROW pointer into Excel output cell array
   header2 {t,1} = upper(m(p).name);
   header2 {t,2} = fs(p);
   header2 {t,3} = output(p).fs;
   header2 {t,4} = lenX(p)/fs(p);
   header2 {t,5} = output(p).lenXtrimmed;
   header2 {t,6} = output(p).activeSpeechLevel;
   header2 {t,7} = output(p).activityFactor;
   header2 {t,8} = output(p).ccMeasures.stDevSum;
   header2 {t,9} = output(p).ccMeasures.varSum;
   header2 {t,10} = output(p).ccMeasures.euclDistMean;
   header2 {t,11} = output(p).ccMeasures.absDistMean;
   header2 {t,12} = output(p).DccMeasures.stDevSum;
   header2 {t,13} = output(p).DccMeasures.varSum;
   header2 {t,14} = output(p).DccMeasures.euclDistMean;
   header2 {t,15} = output(p).DccMeasures.absDistMean;
   header2 {t,16} = output(p).DDccMeasures.stDevSum;
   header2 {t,17} = output(p).DDccMeasures.varSum;
   header2 {t,18} = output(p).DDccMeasures.euclDistMean;
   header2 {t,19} = output(p).DDccMeasures.absDistMean;

end;

% Write output:
xlswrite(outputFile,header1,'Sheet1','H1');
xlswrite(outputFile,header2,'Sheet1','A2');

% Write output to Excel:
header3 = {'Filterbank',' ',' ',' ','Delta Filterbank',' ',' ',' ','Delta Delta Filterbank'};
header4 = {'File','Original fs, Hz','Resampled fs, Hz','Original duration, s',...
   'Trimmed duration, s','Active speech level, dB DFS','Activity factor','sum SD','sum var.','eucl. dist.','abs. dist.'...
   'sum SD','sum var.','eucl. dist.','abs. dist.','sum SD','sum var.','eucl. dist.','abs. dist.'};
for p=1:length(m),
   t = p+1; % ROW pointer into Excel output cell array
   header4 {t,1} = upper(m(p).name);
   header4 {t,2} = fs(p);
   header4 {t,3} = output(p).fs;
   header4 {t,4} = lenX(p)/fs(p);
   header4 {t,5} = output(p).lenXtrimmed;
   header4 {t,6} = output(p).activeSpeechLevel;
   header4 {t,7} = output(p).activityFactor;
   header4 {t,8} = output(p).FBMeasures.stDevSum;
   header4 {t,9} = output(p).FBMeasures.varSum;
   header4 {t,10} = output(p).FBMeasures.euclDistMean;
   header4 {t,11} = output(p).FBMeasures.absDistMean;
   header4 {t,12} = output(p).DFBMeasures.stDevSum;
   header4 {t,13} = output(p).DFBMeasures.varSum;
   header4 {t,14} = output(p).DFBMeasures.euclDistMean;
   header4 {t,15} = output(p).DFBMeasures.absDistMean;
   header4 {t,16} = output(p).DDFBMeasures.stDevSum;
   header4 {t,17} = output(p).DDFBMeasures.varSum;
   header4 {t,18} = output(p).DDFBMeasures.euclDistMean;
   header4 {t,19} = output(p).DDFBMeasures.absDistMean;

end;

% Write output:
xlswrite(outputFile,header3,'Sheet2','H1');
xlswrite(outputFile,header4,'Sheet2','A2');

return;

% Bye!
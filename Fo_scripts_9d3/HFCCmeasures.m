function HFCCmeasures()
% This function prompts for a file or set of files in a directory to process HFCC measures, active
% speech level, and SWIPE pitch statistics. The results are written to an Excel file in the WAV
% directory.

% Mark Skowronski, September 17, 2013

% Prompt for HFCC, SWIPE parameters:
parameters = promptHFCCParameters;
parameters = promptSWIPEParameters(parameters);

% Prompt for files to process:
[fileName,dirName] = uigetfile('*.wav','Select WAV files to process','multiselect','on');
if dirName(1)==0, % cancel button
   return;
end;
if ~iscell(fileName), % single entry, string
   fileName = {fileName}; % convert to cell array for compatibility with multi-file selection
end;

% Init variables:
N = length(fileName);
fs = zeros(1,N); % Hz, sampling rate of each file
lenX = zeros(1,N); % samples, length of untrimmed data in file

% Process each WAV file:
for p=1:N,
   % Read WAV file:
   try
      [x,fs(p)] = wavread([dirName,fileName{p}]);
   catch
      x = randn(5000,1);
      fs(p) = 16000;
   end;
   
   % Store file length:
   lenX(p) = length(x); % samples
   
   % Trim endpoints, calculate HFCC measures:
   [output(p),parameters] = processSentence(x,fs(p),parameters);
end;

% Construct output cell array to sheet 1, standard file info and CC and FB measures:
sheet1 = {'File','Original fs, Hz','Resampled fs, Hz','Original duration, s',...
   'Trimmed duration, s','CC SDS','Delta CC SDS','Delta Delta CC SDS',...
   'FB SDS','Delta FB SDS','Delta Delta FB SDS'};
for p=1:N,
   t = p+1; % ROW pointer into Excel output cell array
   sheet1{t,1} = upper(fileName{p});
   sheet1{t,2} = fs(p);
   sheet1{t,3} = output(p).fs;
   sheet1{t,4} = lenX(p)/fs(p); % sec
   sheet1{t,5} = output(p).lenXtrimmed; % sec
   sheet1{t,6} = output(p).ccMeasures.stDevSum;
   sheet1{t,7} = output(p).DccMeasures.stDevSum;
   sheet1{t,8} = output(p).DDccMeasures.stDevSum;
   sheet1{t,9} = output(p).FBMeasures.stDevSum;
   sheet1{t,10} = output(p).DFBMeasures.stDevSum;
   sheet1{t,11} = output(p).DDFBMeasures.stDevSum;
end;

% Include optional measures:
if parameters.calcASL,
   c = size(sheet1,2); % number of existing columns in sheet 1 output
   header = {'ASL, dB DFS','Activity factor'};
   sheet1(1,c+[1:length(header)]) = header;
   for p=1:N,
      t = p+1;
      sheet1{t,c+1} = output(p).activeSpeechLevel;
      sheet1{t,c+2} = output(p).activityFactor;
   end;
end;
if parameters.calcSWIPE,
   c = size(sheet1,2); % number of existing columns in sheet 1 output
   header = {'P0 mean, Hz','P0 std, Hz','P0 min, Hz','P0 median, Hz','P0 max, Hz',...
             'P0 mean, ST','P0 std, ST','P0 min, ST','P0 median, ST','P0 max, ST','P0 N',...
             'PS mean','PS std','PS min','PS median','PS max'};
   sheet1(1,c+[1:length(header)]) = header;
   for p=1:N,
      t = p+1;
      sheet1{t,c+1} = output(p).SWIPE.statsHz.mean;
      sheet1{t,c+2} = output(p).SWIPE.statsHz.std;
      sheet1{t,c+3} = output(p).SWIPE.statsHz.min;
      sheet1{t,c+4} = output(p).SWIPE.statsHz.median;
      sheet1{t,c+5} = output(p).SWIPE.statsHz.max;
      sheet1{t,c+6} = output(p).SWIPE.statsST.mean;
      sheet1{t,c+7} = output(p).SWIPE.statsST.std;
      sheet1{t,c+8} = output(p).SWIPE.statsST.min;
      sheet1{t,c+9} = output(p).SWIPE.statsST.median;
      sheet1{t,c+10} = output(p).SWIPE.statsST.max;
      sheet1{t,c+11} = output(p).SWIPE.statsHz.N; % same as .statsST.N
      sheet1{t,c+12} = output(p).SWIPE.statsPS.mean;
      sheet1{t,c+13} = output(p).SWIPE.statsPS.std;
      sheet1{t,c+14} = output(p).SWIPE.statsPS.min;
      sheet1{t,c+15} = output(p).SWIPE.statsPS.median;
      sheet1{t,c+16} = output(p).SWIPE.statsPS.max;
   end;
end;

% Delete existing file if necessary:
outputFile = [dirName,'HFCCmeasures.xlsx'];
if exist(outputFile,'file'),
   delete(outputFile);
end;

% Write output:
xlswrite(outputFile,sheet1,'Sheet1','A1');

% Write parameters to Excel file:
sheet2 = {'CC method','Window type','Window size, ms','HFCC Frame rate, fps','FFT size, samples',...
   'Number filters','Number coefficients','Minimum frequency range, Hz','Maximum frequency range, Hz',...
   'E factor','Delta size, frames','Apply pre-emphasis','Pre-emphasis alpha',...
   'Apply cepstral mean subtraction','Resampled sampling rate, Hz','HPF order','HPF cutoff frequency, Hz',...
   'Trim threshold, dB','Trim time constant, ms','Calculate active speech level','Calculate SWIPE',...
   'Minimum pitch range, Hz','Maximum pitch range, Hz','SWIPE Frame rate, fps','P0 candidate spacing, oct',...
   'Spectrum spacing, ERBrate','Window overlap','Pitch strength threshold'};
sheet2(2,:) = {parameters.ccMethod,parameters.windowType,parameters.windowSize*1e3,...
   parameters.frameRate,parameters.FFTsize,parameters.numFilters,parameters.numCoeffs,...
   parameters.freqRange(1),parameters.freqRange(2),parameters.Efactor,parameters.deltaSize,...
   parameters.applyPreemphasis,parameters.preemphasisAlpha,parameters.applyCMS,parameters.fsHFCC,...
   parameters.HPF_order,parameters.HPF_fc,parameters.trimThreshold,parameters.trimTimeConstant*1e3,...
   parameters.calcASL,parameters.calcSWIPE,parameters.plim(1),parameters.plim(2),1/parameters.dt,...
   parameters.dlog2p,parameters.dERBs,parameters.woverlap,parameters.sTHR};

% Write to output:
xlswrite(outputFile,sheet2,'Sheet2','A1');

return;

% Bye!
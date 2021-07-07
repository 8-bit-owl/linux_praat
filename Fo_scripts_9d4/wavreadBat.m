function values = wavreadBat(fileName,N)
% function values = wavreadBat(fileName,N)
% This function reads common .wav or .aiff file formats used in bat or bird experiments.
% Input:  fileName--character array of name of one .wav or .aiff file, including file extension
%         N--samples to input.  If N is scalar, the first N samples are read; if
%            N is a vector, the samples from N(1) to N(2) are read.
% Output: values.x--Nx1 vector of audio data (or Nx2 for stereo, Nx4 for 4-ch data).
%                   8-bit: -128..127, returned as int8 for memory limitations
%                   16-bit: -32768..32767, returned as int16.
%         values.fs--sampling rate, Hz (corrects for time-expanded sampling rates for .wav files)
%         values.source--Estimated source of bat recording (Pettersson, Avisoft, Titley, Sonobat, Batscan, Mac, or unknown)
%         values.bitsPerSample--number of bits per audio sample
%         values.noiseFloordB--estimate of noise floor for int8 or int16 integer data, in dB
%         values.peakSNRdB--estimate of signal-to-noise ratio, in dB: peak energy - noiseFloordB
%         values.audioFormat--1==PCM, 2==Microsoft ADPCM.  Only PCM is currently supported (cannot decode Batscan files).
%                             See http://www.iana.org/assignments/wave-avi-codec-registry for complete list.
%         values.numChannels--number of audio streams in data (1,2, or 4 only)
%         values.riff, .wave, .fmt, .data--sanity check of text present in all .wav headers.  If
%                                          text doesn't match field name, something is wrong...
%         values.lastRIFF--byte number of 'R' in 'RIFF' in header.  Some sources (e.g., Pettersson) have 3 or 4 'RIFF's present.
%         values.chunkSize--total file length (starting from 'R' in last 'RIFF') minus 8 bytes
%         values.subchunk1Size--Number of bytes between 'fmt ' and 'data' text tags (usually 16)
%         values.subchunk2Size--Number of bytes in audio data
%         values.dataPtr: Byte number of 'd' in 'data' tag (offset from lastRIFF value)
% For details on .wav header, see the following: http://ccrma.stanford.edu/CCRMA/Courses/422/projects/WaveFormat/
% Note: values.x = [] for Batscan files because the MS ADPCM codec is not currently supported.

% Version 3, Mark Skowronski, May 7, 2008
% Based on version 2, this file reads .raw files from the batcorder system
% (http://www.ecoobs.de/en/cnt-batcorder.html).

% Version 2, Mark Skowronski, November 28, 2006
% Based on version 1, this file includes support for .aiff audio files from Macs.

% Mark Skowronski, July 6, 2005
% Comments or bug fixes?  Send to mskowro2@uwo.ca

% Check input:
if nargin<2,
   N = [];
end;

% Get file type:
fileExt = lower(fileName(end-2:end));

if strcmp(fileExt,'wav'),
   % Check file for validity:
   fid = fopen(fileName);
   if fid==-1,
      warning(['ERROR: ',fileName,' is not a valid file name.']);
      values(1).lastRIFF = [];
      values.fs = 250e3; % Hz
      values.x = randn(values.fs,1); % simple random noise instead
      values.noiseFloordB = []; % No noise floor estimate
      values.peakSNRdB = []; % No peak SNR estimate
      values.numChannels = 1;
      values.source = 'Noise';
      values.bitsPerSample = 64;
      values.audioFormat = 1;
      return;
   end;

   % Read header:
   fseek(fid,0,-1); % Move cursor to beginning of file
   head = fread(fid,60000,'uint8'); % Read as unsigned bytes
   
   % Init output struct:
   values = struct([]);

   % Find last RIFF:
   h = char(head(:)'); % Convert to ASCII text
   r = findstr('RIFF',h); % At least one, always
   if isempty(r),
      warning(['ERROR: ',fileName,' is not a valid WAV file, no RIFF string in header.']);
      values(1).lastRIFF = [];
      values.fs = 250e3; % Hz
      values.x = randn(values.fs,1); % simple random noise instead
      values.noiseFloordB = []; % No noise floor estimate
      values.peakSNRdB = []; % No peak SNR estimate
      values.numChannels = 1;
      values.source = 'Noise';
      values.bitsPerSample = 64;
      values.audioFormat = 1;
      fclose(fid);
      return;
   end;
   r = r(end); % Point to last one (normal .wav files and Avisoft only have one, but Pettersson have 3 or 4)

   % Point at last RIFF, header follows:
   head = head(r:end);

   % Pull out various parameters:
   values(1).lastRIFF = r;
   values.riff = char(head(1:4))'; % 'RIFF' in ASCII, always
   values.chunkSize = head([1:4]+4); % Chunk size, in little endian
   values.chunkSize = 16.^[0:2:6]*values.chunkSize(:); % Convert 4-byte little endian to number
   values.wave = char(head([1:4]+8))'; % 'WAVE' in ASCII
   values.fmt = char(head([1:4]+12))'; % 'fmt ' in ASCII
   values.subchunk1Size = head([1:4]+16); % 16 bytes for PCM, 18 bytes for Pettersson (2 extra after BitsPerSample), sometimes larger when extra parameters are defined between bitsPerSample and subChunk2ID
   values.subchunk1Size = 16.^[0:2:6]*values.subchunk1Size(:); % Convert 4-byte little endian to number
   values.audioFormat = head([1:2]+20);
   values.audioFormat = 16.^[0,2]*values.audioFormat(:); % Convert 2-byte little endian to number
   values.numChannels = head([1:2]+22);
   values.numChannels = 16.^[0,2]*values.numChannels(:);
   values.fs = head([1:4]+24);
   values.fs = 16.^[0:2:6]*values.fs(:);
%   if values.fs<100e3, % Some bat recording software divides sampling rate by 10
%      values.fs = values.fs*10; % Time-expanded, so correct
%   end;
   values.bitsPerSample = head([1:2]+34);
   values.bitsPerSample = 16.^[0,2]*values.bitsPerSample(:);

   % Find 'data' field (varies since there could be extra parameters present after bitsPerSample):
   dataPtr = findstr('data',char(head(:)')); % one, always (unless also used in a comment field)
   dataPtr = dataPtr(end); % Use last one (could be 'data' in some comments before 'data' field
   values.dataPtr = dataPtr;
   values.data = char(head([1:4]+values.dataPtr-1))';
   values.subchunk2Size = head([1:4]+values.dataPtr-1+4);
   values.subchunk2Size = 16.^[0:2:6]*values.subchunk2Size(:); % Convert 4-byte little endian to number

   if values.audioFormat==1, % PCM (no codec)
      % Get data (8- and 16-bit data treated differently, according to RIFF specification):
      fseek(fid,values.dataPtr-1+values.lastRIFF-1+8,-1); % Move cursor to beginning of data, past 'data' and 4-byte subchunk2
      if values.bitsPerSample==8, % Read as 8-bit unsigned integers, store as int8 variable
         if isempty(N), % Read everything:
            values.x = fread(fid,values.subchunk2Size,'*uint8');
         elseif length(N)==1, % Read samples 1..N/channel:
            if N*values.numChannels<values.subchunk2Size,
               values.x = fread(fid,N*values.numChannels,'*uint8');
            else
               values.x = fread(fid,values.subchunk2Size,'*uint8');
            end;
         else % Read samples N(1)...N(2):
            fseek(fid,values.dataPtr-1+values.lastRIFF-1+8+(N(1)-1)*values.numChannels,-1); % Move cursor to beginning of data
            values.x = fread(fid,(N(2)-N(1)+1)*values.numChannels,'*uint8');
         end;

         % Convert from uint8 (0..255) to int16 (-32768..32767), subtract 128
         % and store as int8 (-128..127).  Note: conversion to int16 is required
         % so that subtracting 128 from the original uint8 variable doesn't saturate.
         values.x = int8(int16(values.x)-128);
      elseif values.bitsPerSample==16, % Read as 16-bit signed integers, return as int16 variable
         if isempty(N), % Read everything:
            values.x = fread(fid,values.subchunk2Size/2,'*int16'); % 16-bit signed integers, 2's complement
         elseif length(N)==1, % Read samples 1..N per channel:
            if N*values.numChannels<values.subchunk2Size/2,
               values.x = fread(fid,N*values.numChannels,'*int16');
            else
               values.x = fread(fid,values.subchunk2Size/2,'*int16');
            end;
         else % Read samples N(1)..N(2):
            fseek(fid,values.dataPtr-1+values.lastRIFF-1+8+2*(N(1)-1)*values.numChannels,-1); % Move cursor to beginning of data
            values.x = fread(fid,(N(2)-N(1)+1)*values.numChannels,'*int16');
         end;
      end;

      % Write interlaced data as multiple channels, if needed:
      if values.numChannels==2,
         values.x = [values.x(1:2:end),values.x(2:2:end)]; % Left, right channels
      elseif values.numChannels==4,
         values.x = [values.x(1:4:end),values.x(2:4:end),values.x(3:4:end),values.x(4:4:end)];
      end;

      % Estimate noise floor:
      if 0,
         fs = values.fs;
         numFrames = floor(length(values.x)/fs/(2/1000)); % 2 ms windows to estimate power
         g = zeros(1,numFrames);
         for p=1:numFrames,
            xIndex = [1:floor(fs*2/1000)]+(p-1)*floor(fs*2/1000); % 2 ms frames, no overlap
            g(p) = mean(double(values.x(xIndex)).^2); % mean squared value
         end;
         g = 10*log10(g+eps); % rms, in dB (eps added in case g(p)=0 for some p)
         [aa,bb] = hist(g,ceil(numFrames/10));
         [junk,bbIndex] = max(aa);
      end;
      values.noiseFloordB = []; % bb(bbIndex); % dB estimate
      values.peakSNRdB = []; % max(g)-bb(bbIndex);
   else % No codec support at this time
      values.x = []; % No data stored
      values.noiseFloordB = []; % No noise floor estimate
      values.peakSNRdB = []; % No peak SNR estimate
   end;

   % Determine source:
   fseek(fid,-100,1); % Move cursor 100 bytes from end of file
   foot = fread(fid,inf,'uint8'); % Read as unsigned bytes
   if ~isempty(findstr('Pettersson',h(:)')), % Pettersson
      values.source = 'Pettersson';
      if values.fs<100e3,
         values.fs = values.fs*10; % Pettersson divides by 10 (sometimes...)
      end;
   elseif ~isempty(findstr('MMMMMMMMM',h(:)')), % Sonobat
      values.source = 'Sonobat';
   elseif ~isempty(findstr('fact',h(:)')) & values.audioFormat==1, % Titley
      values.source = 'Titley'; % 'fact' is normally found when audioFormat>1
   elseif ~isempty(findstr('TIME',foot(:)')) | values.numChannels==4,
      values.source = 'Avisoft';
   elseif values.audioFormat==2, % ADPCM codec, probably from Batbox/Batscan
      values.source = 'Batscan';
   else
      values.source = 'Unknown source';
   end;

   % Close:
   fclose(fid);
elseif strcmp(fileExt,'aif') || strcmp(fileExt,'iff'), % .aif or .aiff file
   % Check file for validity:
   fid = fopen(fileName,'r','b'); % open for reading only, big-endian (Motorola) byte order
   if fid==-1,
      error('ERROR: invalid file name.  Quitting.');
      values(1).lastRIFF = [];
      values.fs = 250e3; % Hz
      values.x = randn(values.fs,1); % simple random noise instead
      values.noiseFloordB = []; % No noise floor estimate
      values.peakSNRdB = []; % No peak SNR estimate
      values.numChannels = 1;
      values.source = 'Noise';
      values.bitsPerSample = 64;
      values.audioFormat = 1;
      return;
   end;

   % Read header:
   fseek(fid,0,-1); % Move cursor to beginning of file
   head = fread(fid,100,'uint8'); % Read as unsigned bytes

   % Find AIFF:
   h = char(head(:)'); % Convert to ASCII text
   if ~strcmp(h(1:4),'FORM') | ~strcmp(h(9:12),'AIFF'),
      error('ERROR: not a valid AIFF file.  Quitting.');
      return;
   end;

   % Init output:
   values = struct([]);

   % Find COMM chunk:
   c = findstr('COMM',h);
   fseek(fid,c+7,-1); % Move file-reading cursor past chunkID and chunkSize
   values(1).numChannels = fread(fid,1,'uint16');
   values.numFrames = fread(fid,1,'uint32');
   values.numBits = fread(fid,1,'uint16');
   values.bitsPerSample = values.numBits;
   ieeeExp = fread(fid,1,'uint16'); % 1-bit sign and 15-bit exponent of 80-bit IEEE float
   ieeeExp = ieeeExp-(2^(15-1)-1); % remove exponent bias
   ieeeMan = fread(fid,1,'uint16'); % most significant 16-bits
   values.fs = ieeeMan*2^(ieeeExp-floor(log(ieeeMan)/log(2))); % Convert mantissa and exponent into double
   values.source = 'Mac';

   % Find SSND chunk:
   c = findstr('SSND',h);
   fseek(fid,c+15,-1); % Move file-reading cursor past chunkID, chunkSize, offset, and blocksize
   if values.numBits==8,
      values.x = fread(fid,inf,'*int8'); % Read as int8, return as int8
   elseif values.numBits==16,
      values.x = fread(fid,inf,'*int16');
   elseif values.numBits==32,
      values.x = fread(fid,inf,'*int32');
   else
      values.x = [];
   end;
   fclose(fid);

   % Estimate noise floor:
   if 0,
      fs = values.fs;
      numFrames = floor(length(values.x)/fs/(2/1000)); % 2 ms windows to estimate power
      g = zeros(1,numFrames);
      for p=1:numFrames,
         xIndex = [1:floor(fs*2/1000)]+(p-1)*floor(fs*2/1000); % 2 ms frames, no overlap
         g(p) = mean(double(values.x(xIndex)).^2); % mean squared value
      end;
      g = 10*log10(g+1e-3); % rms, in dB (eps added in case g(p)=0 for some p)
      [aa,bb] = hist(g,ceil(numFrames/10));
      [junk,bbIndex] = max(aa);
   end;
   values.noiseFloordB = []; % bb(bbIndex); % dB estimate
   values.peakSNRdB = []; % max(g)-bb(bbIndex);
   
   values.audioFormat = 1;
   values.numChannels = 1;
elseif strcmp(fileExt,'raw'), % batcorder file, 500 kHz, 16 bit, little endian
   % Open file:
   fid = fopen(fileName,'r','ieee-le'); % read only, little-endian byte order
   if fid==-1,
      error('ERROR: invalid file name.  Quitting.');
      values(1).lastRIFF = [];
      values.fs = 250e3; % Hz
      values.x = randn(values.fs,1); % simple random noise instead
      values.noiseFloordB = []; % No noise floor estimate
      values.peakSNRdB = []; % No peak SNR estimate
      values.numChannels = 1;
      values.source = 'Noise';
      values.bitsPerSample = 64;
      values.audioFormat = 1;
      return;
   end;
   
   % Init output:
   values = struct([]);

   % Find COMM chunk:
   values(1).numChannels = 1;
   values.numBits = 16;
   values.bitsPerSample = values.numBits;
   values.fs = 500e3; % Hz
   values.source = 'batcorder';
   
   % Read entire file:
   fseek(fid,0,-1); % Move cursor to beginning of file
   values.x = fread(fid,inf,'*int16'); % 16-bit signed, all data, return as 16-bit data
   
   values.noiseFloordB = []; % bb(bbIndex); % dB estimate
   values.peakSNRdB = []; % max(g)-bb(bbIndex);
   
   values.audioFormat = 1;
   values.numChannels = 1;
else
   values = struct([]);
   values(1).lastRIFF = [];
   values.fs = 250e3; % Hz
   values.x = randn(values.fs,1); % simple random noise instead
   values.noiseFloordB = []; % No noise floor estimate
   values.peakSNRdB = []; % No peak SNR estimate
   values.numChannels = 1;
   values.source = 'Noise';
   values.bitsPerSample = 64;
   values.audioFormat = 1;
end;


% Bye!

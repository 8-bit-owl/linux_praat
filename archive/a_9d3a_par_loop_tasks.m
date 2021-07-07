% control_loop7d8h
%
% loads in wave files from a directory to process for dB and Fo and other
% stuff.  it has been optimized for long files from an accelerometer but
% could work from microphone data.  it requires that the following three files
% are in the same directory as this script:
%   praat_pitch.psc
%   praat_voice.psc
%   praatcon.exe (new version is just praat.exe)
%
% two different F0 extraction techniques are built in.
% Also, accompaning scripts are necesar in a path location:
%   R:\Data\DATA-priority03-Research\Tobias-Bird\Accelerometer\Fo_scripts
% These scripts in this folder are required.  the path shoudl be updated to
% match whereever you housetue scripts.
%
% eric.hunter
% 2013-June
% 2013-Oct -  couple small edits, export Fo Mode
% 2014-clear allJan -  edit to open full wave file.
% 2014-Feb -  added pitch strength using AudSwipe
% 2014-Mar 6a, add calibration and dB stats
% 2014-Mar 6b, mark parameters
% 2014-M                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ```````````````````````````````````````````````````````````````````````````````ay 6c2, only do one channel but you can name it, added new swipe.
% 2014-May 6d, added voice measures to the analysis
% 2014-July 6e, reverted to original AudSwipe, and streamlined plots
% 2014-July 6f, adjusted the duration find and used it for the spectral
% analysis.  see Jitter and Shimmer.  Also, added concatenated.
% 2014-AUg 6g, adjusted to better get voicing and speech
% 2014-AUg 6h, bug fix on concat and output
% 2014-Spt 6i, adjusted concat to have an option just to keep the data but
% not be too long by keeping waveforms
% 2014-Oct 6j, adjusted the LTAS because of error.
% 2014-Nov 6j2, slight adjustments to check output
% 2014-Dec 6j3, make to be both task and steady vowel, add level threshold
% 2015-Oct 6j4, update to change in PRAATCON_win98.exe (praat updated away
% from the command line but it didnt work with what we were doing)
% 2016-Feb 6k1, Mark Berardi - multiple changes
% 2016-Nov 7d2, - FilenamesByExt change
% 2018-July 7f2, EJH, updated the dur0task estimation as well as some minor
% cleanups
% 2018-Aug 7f4, ejh, added fry analysis for modal and other numbers, also
% adjusted the way some of the results are Concatenated
% 2018-Aug 7f5, ejh, added generic script for praat analysis, this allows
% for the range of fo to be based solely on fo_lower, fo_upper for both
% praat and audswipe. for praat, a temporary file is written for the
% analysis which also has a unique identifier so that mutliple can use the
% same fo script folder.
% 2018-Aug 7f6, ejh adjust the voice praat so that multiple can get it.
% 2018-sept 7f7, add the ability to export fo and db contours
% 2018-sept 7f7c, updated AVQI to export cpps and related for both full
% sample and the voice only part.
% 2018-sept 7f7d, updated to add voice detection variable to pratt extraction
% 2019-Jan 7f7e, updated plots of audswipe showing 0s
% 2019-Mar 7f7f, bug fixes, updated entropy and file names
% 2019-Apr 7f8, added fry percentage on output
% 2019-Apr 7f8g, update praat scripts to be updated from commandline
% 2019-May 8hb, update output to excel. added name to files, auto Fo range find
% 2019-May 8hc, add SHRP as one more extraction method and other updates
% 2019-May 8hd, added new metrics and also updated output and graphics
% 2019-May 9, stable with new additions
% 2019-May par9, same analysis but with parallel
% 2019-May par9a, update for vowels
% 2019-May par9b, updates some errors
% 2019-May par9c, fixed errors, updated analy folder, updated csv output

% clear all;
% close all;
% clc;
tic

%% settings
   % 0 sets the setting to default
s_dir = 0;   % 0: use default path.
   % 1: use a code location and file location (set in next section)
   % 2: be prompted by a menu to choose file locations
s_batch = 0; % 0: use FilenameByExt.m to batch files
   % 1: use a menu to select files to batch together
s_gen = 0;   % 0: select user pitch rages (must change all code/subroutines).
   % 1: uses ranges for female voices (150-800 Hz)
   % 2: uses ranges for male voices (70-350 Hz)
   % 3: uses ranges for both voices depending on gennum
s_gennum = 0; % 0: no gennum needed
   % #: > 0 number that is the number character in the filename
   % that indicated gender (m or M or f or F)
s_name = 0;  % 0: use default file output name
   % 1: will be prompted to name the output .xls file

s_start = 1;  % if crash, put j_iter+1 here and continue running (F5), note bad file


%% initialize USER EDITABLE
u_ftitle = 'MCDDK'; %'% name to put in documents, no spaces';
u_avqi = 0;  % 0: do not run AVQI analysis
   % 1: run both AVQI analyses, be sure to change the AVQI settings below
u_fftenvelope = 0; % 0 to skip, 1 to do the fft of the array or envelope to find features.
u_trimYN = 0;         % 0 to analyze the files as they are
   % 1 use dB as a segmenter
   % 2 to do speech detection and toss small stuff
   % 3 keep only voicing, good for running speech
   % 4 find beginning and end of voicing/dB and just do that, good for steady
   % vowels that have no other things.
u_concatYN = 1;     % 2 to concatenate and also save the concatenated wave file
   % 1 to concatenate the voicing sections,
   % 0 to only do individual analysis
u_contour2xlsYN = 1;     % 1 to export the arrays to a spreadsheet (Fo, dB etc)
   % 0 does not.

c_fryanaly = 0; % 0 everything is analyzed. 1 only modal voice is analyzed, 2 only fry is analyzed
c_autofo = 0; % 0: means that just do the values set below
   % 1: will do an auto find range for the file
c_foS_low=140;   %lower bound for Fo (70-90 for males, 150 for females)
c_foS_up=370;   %upper bound for Fo (350 for males, 800 for females)
   % fo_lower=35;    %lower bound for Fo (70-90 for males, 150 for females)
   % fo_upper=85;   %upper bound for Fo (350 for males, 800 for females)
u_showfindFofig = 5; % 0 to not show the figure, 5 to show and save
u_figson = 1;       % no monitoring figures =0
   % monitoring each signal being loaded = 1
   % monitoring the output as it is calculated = 2
   % monitoring each signal and output as is calculated = 3

c_LTAS_lower = 50;    % lower range to show LTAS, calculate slopes, etc
c_LTAS_upper = 5000; % upper range to show LTAS, calculate slopes, etc

% the above also need to be adjusted in file "praat_pitch.psc" if changes
% are made.  [To Pitch (ac)... 0 80 15 no 0.03 0.45 0.0025 0.35 0.20 700]
%  [To Pitch (ac)... 0 [fo_lower] 15 no 0.03 0.45 0.0025 0.35 0.20 [fo_upper]]

if s_gen == 1
   c_foS_low=150;    %lower bound for Fo (70-90 for males, 150 for females)
   c_foS_up=800;   %upper bound for Fo (350 for males, 800 for females)
elseif s_gen == 2
   c_foS_low=65;    %lower bound for Fo (70-90 for males, 150 for females)
   c_foS_up=350;   %upper bound for Fo (350 for males, 800 for females)
end
s_g0 = s_gen;
set(0,'DefaultFigureWindowStyle','docked')

s_micOacc=1;  % enter 1 for mic and 2 for acc
s_gain_voltA=20000; % Gain voltage for the recording
% gain_voltA=450025.1; % Gain voltage for the recording


%% Get Path information (sdir)
if s_dir == 0
   %code = 'V:\1 Software - Equipment - Techniques -Analysis\Matlab Universal Scripts\Fo_scripts';    % network
   s_code = 'C:\cloud\OneDrive - Michigan State University\Current Work\1 Matlab analysis\Fo_scripts_9d3';    % local
%    s_code = 'C:\1analyC\Fo_scripts_9d3';  
   s_files = cd;
   addpath(s_code)
   
elseif s_dir == 1
   s_files = 'V:\Offsite Data Collection\LSVT loud\redo low voice';
   s_code = 'V:\1 Software - Equipment - Techniques -Analysis\Matlab Universal Scripts\Fo_scripts';
   addpath(s_code)
   
elseif s_dir == 2
   s_files = uigetdir(cd,'Choose file directory');
   s_code = 'V:\1 Software - Equipment - Techniques -Analysis\Matlab Universal Scripts\Fo_scripts';
   addpath(s_code)
end


%% get file list to analyze (sbatch)
if s_batch == 0
   oldfolder=cd(s_files);
   %     filename = FilenamesByExt('wav',files);
   filetmp1= FilenamesByExt('wav');
   filetmp2 = FilenamesByExt('WAV');
   filetmp3 = FilenamesByExt('Wav');
   g_filename=[filetmp1(:)',filetmp2(:)',filetmp3(:)']';
   cd(oldfolder)
   clear filetmp1 filetmp2 filetmp3 oldfolder
elseif s_batch == 1
   g_filename = uigetfile({'*.wav;*.WAV;*.Wav'},'Choose wav files to analyze',s_files,'MultiSelect','on');
end

g_cnt=length(g_filename);

%% %%%%%%%% dont change below
if s_micOacc==1, g_titlename=['_' u_ftitle '_1MIC']; else, g_titlename=['_' u_ftitle '_1ACC']; end
if c_fryanaly == 1; g_titlename= [g_titlename '-mod']; end
if c_fryanaly == 2; g_titlename= [g_titlename '-fry']; end
% skp=25;           %how many points to skip in plots
s_RMS_window=20;    %msec
% fo_window=20;     %window size for getting Fo (in msec)
s_foS_step=10;       %step of the Fo window (in msec)
% PlayBackSpeed=1.0;

s_foS_Analysis = 'ac'; % 'cc'
s_foS_NCand = 15; % max number of candidates
s_foS_Accuracy = 'no';  % 'no'
s_foS_SilenceThrsh = 0.03; % silenceThresh - (standard value: 0.03); frames that do not contain amplitudes above this threshold (relative to the global maximum amplitude), are probably silent.
s_foS_VoiceThrsh= 0.4; %voicing threshold - (standard value: 0.45); the strength of the unvoiced candidate, relative to the maximum possible autocorrelation. To increase the number of unvoiced decisions, increase this value.
s_foS_OctCost =0.01; % octave cost - (standard value: 0.01 per octave); degree of favouring of high-frequency candidates, relative to the maximum possible autocorrelation. This is necessary because even (or: especially) in the case of a perfectly periodic signal, all undertones of F0 are equally strong candidates as F0 itself. To more strongly favour recruitment of high-frequency candidates, increase this value.
s_foS_OctJumpCost = 0.35; % octave-jump cost - (standard value: 0.35); degree of disfavouring of pitch changes, relative to the maximum possible autocorrelation. To decrease the number of large frequency jumps, increase this value. In contrast with what is described in the article, this value will be corrected for the time step: multiply by 0.01 s / TimeStep to get the value in the way it is used in the formulas in the article.
s_foS_VoiceUnvoiceCost = 0.25; % voiced/unvoiced cost - (standard value: 0.14); degree of disfavouring of voiced/unvoiced transitions, relative to the maximum possible autocorrelation. To decrease the number of voiced/unvoiced transitions, increase this value. In contrast with what is described in the article, this value will be corrected for the time step: multiply by 0.01 s / TimeStep to get the value in the way it is used in the formulas in the article.
s_pertS_MaxPeriodFact=1.3;
s_pertS_MaxAmpFact=1.6;
% description of some of the commands is here: http://www.fon.hum.uva.nl/praat/manual/Voice_2__Jitter.html

%% variables for concatination
% allocate initial variables to save the data in

cncat1_sigS = []; %concatenated audio signal
cncat1_sigS_t = [];
cncat1_sig_dB = [];
cncat1_sig_dBs = [];
cncat1_sig_dBi = [];
cncat1_t_t = [];
cncat1_s_t = [];
cncat1_i_t = [];
% cncat1_t_tc = [];
% cncat1_s_tc = [];
% cncat1_i_tc = [];
cncat1_f0_So_aud = [];
cncat1_f0_So_auds = [];
cncat1_f0_So_audi = [];
cncat1_f0_fo_shp = [];
cncat1_f0_fo_shps = [];
cncat1_f0_fo_shpi = [];
cncat1_f0_fo_prtf = [];
cncat1_f0_fo_prtfs = [];
cncat1_f0_fo_prtfi = [];
cncat1_f0_fo_aud = [];
cncat1_f0_fo_auds = [];
cncat1_f0_fo_audi = [];
cncat1_g_duration = [];

% Output name
if s_name == 0
   if s_micOacc==1, g_titlename=['_' u_ftitle '_1MIC']; else, g_titlename=['_' u_ftitle '_1ACC']; end
   if c_fryanaly == 1; g_titlename= [g_titlename '-mod']; end
   if c_fryanaly == 2; g_titlename= [g_titlename '-fry']; end
   g_outname = ['\',g_titlename,'_FullResults'];
elseif s_name == 1
   tmp = inputdlg('Name of output file (without ''.xls'')');
   g_outname = ['\' tmp{1}];
end
clear tmp

% file writes
fclose('all'); 
if u_contour2xlsYN==1
   g_fid_out = fopen([g_titlename '_contour.csv'],'w'); %open file and create fid
   fprintf(g_fid_out,'%s,\t%s,\t%s,\t%s,\n', 'Filename','Metrics','Subdetails', 'data___');
   fclose(g_fid_out); %close file
end
delete([g_titlename '_trackcomplete.txt']);
g_fid_out = fopen([g_titlename '_trackcomplete.txt'],'w');
fclose(g_fid_out); %close file
delete([g_titlename '_failed.txt']);
g_fid_out = fopen([g_titlename '_failed.txt'],'w');
fclose(g_fid_out); %close file

%% primary loop
% for g_iter = s_start:g_cnt          %loop through all files in the directory
parfor g_iter = s_start:g_cnt          %loop through all files in the directory
   %%
   try
   [outdata,sigaS,t1S,sig_dB,t2,t2Si,dB_vals,f0_So_aud,f0_fo_shp,f0_fo_prtf,f0_fo_aud,...
      g_duration]=...
      a_9d3a_par_task_function(g_iter,g_cnt,g_filename,g_titlename,...
      s_gen,s_micOacc,s_gennum,u_avqi,s_gain_voltA,s_code,s_files,...
      s_RMS_window,s_foS_step,s_foS_Analysis,s_foS_NCand,s_foS_Accuracy,s_foS_SilenceThrsh,...
      s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,s_foS_VoiceUnvoiceCost,...
      s_pertS_MaxPeriodFact,s_pertS_MaxAmpFact,...
      u_trimYN,u_concatYN,u_contour2xlsYN,u_fftenvelope,u_figson,u_showfindFofig,...
      c_LTAS_lower,c_LTAS_upper,c_autofo,c_foS_low,c_foS_up,c_fryanaly);
    outdataN(g_iter,:)=outdata;
   
   % 1_______________ Concatinate _____________
   if u_concatYN>0
      if u_concatYN==2
         cncat1_sigS = [cncat1_sigS sigaS'];
         cncat1_sigS_t = [cncat1_sigS_t t1S];
      end
      cncat1_g_duration=[cncat1_g_duration g_duration];
      tmp=find(t2Si>0);
      cncat1_sig_dB = [cncat1_sig_dB sig_dB];    % concatenate the dB
      cncat1_sig_dBs = [cncat1_sig_dBs sig_dB(tmp(1):tmp(end))];    % concatenate dB of just detected speech
      cncat1_sig_dBi = [cncat1_sig_dBi sig_dB(t2Si)];    % concatenate the dB
      
      cncat1_t_t = [cncat1_t_t t2/60]; % change to minutes
      cncat1_s_t = [cncat1_s_t t2(tmp(1):tmp(end))/60]; % change to minutes
      cncat1_i_t = [cncat1_i_t t2(t2Si)/60]; % change to minutes
      
      %       cncat1_t_tc = [cncat1_t_tc cncat1_t_t(end)+t2/60]; % change to minutes
      %       cncat1_s_tc = [cncat1_s_tc cncat1_s_t(end)+t2(tmp(1):tmp(end))/60]; % change to minutes
      %       cncat1_i_tc = [cncat1_i_tc cncat1_i_t(end)+t2(t2Si)/60]; % change to minutes
      
      cncat1_f0_So_aud = [cncat1_f0_So_aud f0_So_aud'];
      cncat1_f0_So_auds = [cncat1_f0_So_auds f0_So_aud(tmp(1):tmp(end))'];
      cncat1_f0_So_audi = [cncat1_f0_So_audi f0_So_aud(t2Si)'];
      
      cncat1_f0_fo_shp = [cncat1_f0_fo_shp f0_fo_shp'];
      cncat1_f0_fo_shps = [cncat1_f0_fo_shps f0_fo_shp(tmp(1):tmp(end))'];
      cncat1_f0_fo_shpi = [cncat1_f0_fo_shpi f0_fo_shp(t2Si)'];
      
      cncat1_f0_fo_prtf = [cncat1_f0_fo_prtf f0_fo_prtf'];
      cncat1_f0_fo_prtfs = [cncat1_f0_fo_prtfs f0_fo_prtf(tmp(1):tmp(end))'];
      cncat1_f0_fo_prtfi = [cncat1_f0_fo_prtfi f0_fo_prtf(t2Si)'];
      
      cncat1_f0_fo_aud = [cncat1_f0_fo_aud f0_fo_aud'];
      cncat1_f0_fo_auds = [cncat1_f0_fo_auds f0_fo_aud(tmp(1):tmp(end))'];
      cncat1_f0_fo_audi = [cncat1_f0_fo_audi f0_fo_aud(t2Si)'];      
      %       clear tmp
   end
   catch  
      %% write out to success file
      g_fid_out = fopen([g_titlename '_trackcomplete.txt'],'a'); %open file and create fid
      while g_fid_out<0
         fprintf('...waiting...')
         g_fid_out = fopen([g_titlename '_trackcomplete.txt'],'a'); %open file and create fid
      end
      fprintf(g_fid_out,'-%d\n',g_iter); %write out title
      fclose(g_fid_out); %close file
      %% write out failed files
      g_fid_out = fopen([g_titlename '_failed.txt'],'a'); %open file and create fid
      while g_fid_out<0
         fprintf('...waiting...')
         g_fid_out = fopen([g_titlename '_failed.txt'],'a'); %open file and create fid
      end
      fprintf(g_fid_out,'%s\n',g_filename{g_iter}); %write out title
      fclose(g_fid_out); %close file
   end
   
end

fprintf('\n')
%% write output

outdataD='----------';    % metrics names
outdataD=[outdataD;'*SETTINGS*';...
   'trimOption'; 'fry_analy_'; 'RMS_windw_'; 'Fo_WinStp_';...
   'perSilncTh'; 'perVoicTh_'; 'perOctCst_'; 'perOctJmp_';...
   'perVcUnVcC'; 'perMxPerF_'; 'perMxAmpF_'; ...
   'Fo_ExtProb'; 'Fo_LwrBndr'; 'Fo_uprBndr'; 'LTASlwr_Bd'; 'LTASupprBd'];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*GEN_DATA_';...
   'DurBeg2End';...      % end voicing - beginning voicing
   'DurSpchDtc';...      % estimated length of all speaking parts
   'DurVoicALL';...      % duration of voice only (fo segments from all 3 extraction methods)
   'DurVoicPRT';...      % duration of voice only (Praat only voicing)
   '%VoicCnSpc';...      % Voice % compared to the concatenated speech
   '%VoicAlSpc';...      % Voice % from the length of speech, from first instance to last instance
   'FileTime__'];
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*GEN_CALC_';...
   'ActvtyFact';...      % activity factor activityFactor [0 1]
   'FryPerc   ';...      % amount of the recording which is fry compared to the modal segments
   'AlphaRatio'; 'dBspecSlp_'; 'dBspecSlpF'; 'dBspcSpOct';... % alpha ratio, spec slope from 50hz, spec slope from Fo, slop using 3rd oct
   'dB1kTo3k  ';...      % energy from 1kHz to 3kHz normed to that range
   'dB1kTo3kN ';...      % percent energy from 1kHz to 3kHz compared to 50hz-10khz
   'dBSngFmnt ';...      % energy from 2.8kHz to 3.4kHz (singers formant area)
   'dBSngFmntN';...      % percent energy from 2.8kHz to 3.4kHz (singers formant area) to LTAS area
   'HFCCSDS   ';...      % HFCC HFCC standard deviation sum measure
   'DHFCCSDS  ';...      % delta HFCC features, standard deviation sum measure
   'DDHFCSDS  '];        % delta delta HFCC features, standard deviation sum measure
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*CPPctSPCH';...
   'cpp__avCAT';'cpps_avCAT';'hnr_avCAT ';'shim_avCAT';'shdb_avCAT';...
   'slop_avCAT';'tilt_avCAT';'avqi_avCAT'];
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*CPPallFIL';...
   'cpp__avALL';'cpps_avALL';'hnr_avALL ';'shim_avALL';'shdb_avALL';...
   'slop_avALL';'tilt_avALL';'avqi_avALL'];
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*CALC_PURT';...
   'jitter    ';'jitt_abs  ';'jitt_rap  ';'jit_ppq5  ';'jit_ddp   ';...
   'shimmer   ';'shim_db   ';'shim_apq3 ';'shim_apq5 ';...
   'shim_apq11';'shim_dda  ';...
   'nhr       ';'hnr       ';...
   'rpde      ';'dfa       ';'ppe       ';...
   'purt_dur  '];
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*dB_ALLraw';...
   'dB__N     '; 'dB__mean  '; 'dB__std   '; 'dB__var   ';...
   'dB__mode  '; 'dB__med   '; 'dB__Q1    '; 'dB__Q3    ';...
   'dB__IQR   '; 'dB__skew  '; 'dB__bowskw'; 'dB__kurt  ';...
   'dB__slop  '; 'dB__fftslp'; 'dB__ffthz1'; 'dB__ffthz2'; 'dB__ffthz3'; 'dB__ffthz4';...
   'dB__fftA1 '; 'dB__fftA2 '; 'dB__fftA3 '; 'dB__fftA4 '];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*dBvALLspc';...
   'dBv_N     '; 'dBv_mean  '; 'dBv_std   '; 'dBv_var   ';...
   'dBv_mode  '; 'dBv_med   '; 'dBv_Q1    '; 'dBv_Q3    ';...
   'dBv_IQR   '; 'dBv_skew  '; 'dBv_bowskw'; 'dBv_kurt  ';...
   'dBv_slop  '; 'dBv_fftslp'; 'dBv_ffthz1'; 'dBv_ffthz2'; 'dBv_ffthz3'; 'dBv_ffthz4';...
   'dBv_Afft1 '; 'dBv_Afft2 '; 'dBv_Afft3 '; 'dBv_Afft4 '];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Sov_PchSt';...
   'Sov_N     '; 'Sov_mean  '; 'Sov_std   '; 'Sov_var   ';...
   'Sov_mode  '; 'Sov_med   '; 'Sov_Q1    '; 'Sov_Q3    ';...
   'Sov_IQR   '; 'Sov_skew  '; 'Sov_bowskw'; 'Sov_kurt  ';...
   'Sov_slop  '; 'Sov_fftslp'; 'Sov_ffthz1'; 'Sov_ffthz2'; 'Sov_ffthz3'; 'Sov_ffthz4';...
   'Sov_Afft1 '; 'Sov_Afft2 '; 'Sov_Afft3 '; 'Sov_Afft4 '];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Fo_Hz_ALL';...
   'FoAl_N    '; 'FoAl_mean '; 'FoAl_std  '; 'FoAl_var  ';...
   'FoAl_mode '; 'FoAl_med  '; 'FoAl_Q1   '; 'FoAl_Q3   ';...
   'FoAl_IQR  '; 'FoAl_skew '; 'FoAlbowskw'; 'FoAl_kurt ';...
   'FoAl_slop '; 'FoAlfftslp'; 'FoAlffthz1'; 'FoAlffthz2'; 'FoAlffthz3'; 'FoAlffthz4';...
   'FoAl_Afft1'; 'FoAl_Afft2'; 'FoAl_Afft3'; 'FoAl_Afft4';...
   'FoAlmeanST'; 'FoAmnpsdST'; 'FoAmnmstST'];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Fo_ST_ALL';...
   'FoAl_N    '; 'FoAl_mean '; 'FoAl_std  '; 'FoAl_var  ';...
   'FoAl_mode '; 'FoAl_med  '; 'FoAl_Q1   '; 'FoAl_Q3   ';...
   'FoAl_IQR  '; 'FoAl_skew '; 'FoAlbowskw'; 'FoAl_kurt ';...
   'FoAl_slop '];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Fov_Hz___';...
   'Fov_N     '; 'Fov_mean  '; 'Fov_std   '; 'Fov_var   ';...
   'Fov_mode  '; 'Fov_med   '; 'Fov_Q1    '; 'Fov_Q3    ';...
   'Fov_IQR   '; 'Fov_skew  '; 'Fov_bowskw'; 'Fov_kurt  ';...
   'Fov_slop  '; 'Fov_fftslp'; 'Fov_ffthz1'; 'Fov_ffthz2'; 'Fov_ffthz3'; 'Fov_ffthz4';...
   'Fov_Afft1 '; 'Fov_Afft2 '; 'Fov_Afft3 '; 'Fov_Afft4 ';...
   'Fov_meanST'; 'FovmnpsdST'; 'FovmnmsdST'];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Fov_ST___';...
   'Fov_N     '; 'Fov_mean  '; 'Fov_std   '; 'Fov_var   ';...
   'Fov_mode  '; 'Fov_med   '; 'Fov_Q1    '; 'Fov_Q3    ';...
   'Fov_IQR   '; 'Fov_skew  '; 'Fov_bowskw'; 'Fov_kurt  ';...
   'Fov_slop  '];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Pov_Hz___';...
   'Pov_N     '; 'Pov_mean  '; 'Pov_std   '; 'Pov_var   ';...
   'Pov_mode  '; 'Pov_med   '; 'Pov_Q1    '; 'Pov_Q3    ';...
   'Pov_IQR   '; 'Pov_skew  '; 'Pov_bowskw'; 'Pov_kurt  ';...
   'Pov_slop  '; 'Pov_fftslp'; 'Pov_ffthz1'; 'Pov_ffthz2'; 'Pov_ffthz3'; 'Pov_ffthz4';...
   'Pov_Afft1 '; 'Pov_Afft2 '; 'Pov_Afft3 '; 'Pov_Afft4 ';...
   'Pov_meanST'; 'PovmnpsdST'; 'PovmnmsdST'];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Pov_ST___';...
   'Pov_N     '; 'Pov_mean  '; 'Pov_std   '; 'Pov_var   ';...
   'Pov_mode  '; 'Pov_med   '; 'Pov_Q1    '; 'Pov_Q3    ';...
   'Pov_IQR   '; 'Pov_skew  '; 'Pov_bowskw'; 'Pov_kurt  ';...
   'Pov_slop  '];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Fo_SHR_Hz';...
   'FoS_N     '; 'FoS_mean  '; 'FoS_std   '; 'FoS_var   ';...
   'FoS_mode  '; 'FoS_med   '; 'FoS_Q1    '; 'FoS_Q3    ';...
   'FoS_IQR   '; 'FoS_skew  '; 'FoS_bowskw'; 'FoS_kurt  ';...
   'FoS_slop  '; 'FoS_fftslp'; 'FoS_ffthz1'; 'FoS_ffthz2'; 'FoS_ffthz3'; 'FoS_ffthz4';...
   'FoS_Afft1 '; 'FoS_Afft2 '; 'FoS_Afft3 '; 'FoS_Afft4 ';...
   'FoS_meanST'; 'FoSmnpsdST'; 'FoSmnmsdST'];%slope accross utterance
outdataD=[outdataD;'----------'];
outdataD=[outdataD;'*Fo_SHR_ST';...
   'FoS_N     '; 'FoS_mean  '; 'FoS_std   '; 'FoS_var   ';...
   'FoS_mode  '; 'FoS_med   '; 'FoS_Q1    '; 'FoS_Q3    ';...
   'FoS_IQR   '; 'FoS_skew  '; 'FoS_bowskw'; 'FoS_kurt  ';...
   'FoS_slop  '];%slope accross utterance
outdataD=[outdataD;'----------'];


disp(' --start-- Save Gen Results')
outdataN=double(outdataN);
save([s_files,g_outname,'_mat.mat'],'g_cnt','g_filename','outdataD','outdataN')
save([g_titlename '_FullResults.txt'],'outdataN','-ascii')
disp(' ')
disp(' --start-- Listing Output Parameters')
disp(outdataD)
disp(' --end-- Listing Output Parameters')
pause(0.5)
disp(' ')
disp(' --start-- Listing filenames in order of analysis for Excel')
disp(g_filename)
disp(' ')
disp(' --end-- Listing filenames in order of analysis for Excel')

%% Save output to file:
% workSheet = 'datanames';
% g_outname=[g_outname '_' titlename];
workSheet = 'data';
xlswrite([s_files,g_outname,'.xlsx'],[{'Filenames'} cellstr(outdataD)'],workSheet,'A1');
xlswrite([s_files,g_outname,'.xlsx'],outdataN,workSheet,'B2');
xlswrite([s_files,g_outname,'.xlsx'],g_filename,workSheet,'A2');



%% Concatenated Analysis

if u_concatYN>0
   
   disp(' --start-- Concatenated analysis')
   %% Concatinated results
   concat2=[];
   concat2.files=g_filename;
   concat2.g_durationavg=mean(cncat1_g_duration);
   concat2.g_durationtot=sum(cncat1_g_duration);
   
   %     tmp=cncat1_f0_value1_t(~isnan(cncat1_f0_value1));
   rngtmp=s_foS_step/1000;
   concat2.durOspch1=length(cncat1_s_t)*s_foS_step/1000;
   concat2.durOspch2=length(cncat1_i_t)*s_foS_step/1000;
   tmp=[cncat1_f0_fo_shpi cncat1_f0_fo_prtfi cncat1_f0_fo_audi];
   concat2.durOvoice= s_foS_step*sum(~isnan(tmp))/3/1000;
   concat2.durOvoicePr= s_foS_step*sum(~isnan(cncat1_f0_fo_prtf))/1000;
   concat2.perOvoice1= 100*concat2.durOvoice/concat2.durOspch2;
   concat2.perOvoice2= 100*concat2.durOvoice/concat2.durOspch1;
   
   [concat2.FoAllhz_stats,concat2.FoAllst_stats] = ...
      freq_stats([cncat1_f0_fo_prtfi cncat1_f0_fo_shpi cncat1_f0_fo_audi],...
      [cncat1_i_t cncat1_i_t cncat1_i_t],220);
   [concat2.FoPRThz_stats,concat2.FoPRTst_stats] = freq_stats(cncat1_f0_fo_prtfi',cncat1_i_t',220);
   [concat2.FoSHRhz_stats,concat2.FoSHRst_stats] = freq_stats(cncat1_f0_fo_shpi',cncat1_i_t',220);
   [concat2.FoAUDhz_stats,concat2.FoAUDst_stats] = freq_stats(cncat1_f0_fo_audi',cncat1_i_t',220);
   
   tmp=cncat1_f0_So_audi;tmp=tmp(~isnan(tmp)); % AudSwipe Ps
   tmpt=cncat1_i_t; tmpt=tmpt(~isnan(tmp)); p = polyfit(tmpt,tmp,1);
   concat2.SoAud_stats= basicstats(round(tmp)); concat2.SoAud_stats.slope=p;
   
   concat2.dB_stats= basicstats(cncat1_sig_dBs);
   p = polyfit(cncat1_s_t,cncat1_sig_dBs,1);  concat2.dB_stats.slope=p;
   
   concat2.dBv_stats= basicstats(cncat1_sig_dBi);
   p = polyfit(cncat1_i_t,cncat1_sig_dBi,1); concat2.dBv_stats.slope=p;
   
   %    if u_concatYN==2
   %       concat2.jitter1 = voice_measures2(cncat1_sigS, g_Fs);
   %       concat2.jitter1.length=length(cncat1_sigS_t)*(t1(2)-t1(1));
   %       HFCCtmp.x=cncat1_sigS;
   %       [HFCCtmp,HFCCparm] = processSentences(HFCCtmp,HFCCparm);
   %       concat2.c_HFCC.HFCC_SDS=HFCCtmp.ccMeasures.stDevSum;
   %       concat2.c_HFCC.DHFCC_SDS=HFCCtmp.DccMeasures.stDevSum;
   %       concat2.c_HFCC.DDHFCC_SDS=HFCCtmp.DDccMeasures.stDevSum;
   %       concat2.LTAS=ltas(cncat1_sigS,g_Fs,min([2^floor(log2(length(cncat1_sigS))) 2^13]),'han',0,0,1,-100);
   %       f1=concat2.LTAS.f;
   %       spec1=concat2.LTAS.linspectrum;
   %       spec1dB=concat2.LTAS.norm0spectrum;
   %
   %       concat2.LTAS.AlphaRatio=...
   %          sum(concat2.LTAS.norm0spectrum((concat2.LTAS.f>1000&concat2.LTAS.f<=5000)))/sum((concat2.LTAS.f>1000&concat2.LTAS.f<5000))-...
   %          sum(concat2.LTAS.norm0spectrum((concat2.LTAS.f>50&concat2.LTAS.f<=1000)))/sum((concat2.LTAS.f>50&concat2.LTAS.f<1000));
   %       concat2.LTAS.dB1kTo3k=...
   %          sum(concat2.LTAS.norm0spectrum((concat2.LTAS.f>=1000&concat2.LTAS.f<=3000)))...
   %          /sum(concat2.LTAS.f>=1000&concat2.LTAS.f<=3000);
   %       concat2.LTAS.dB1kTo3kN=...
   %          sum(concat2.LTAS.norm0spectrum((concat2.LTAS.f>=1000&concat2.LTAS.f<=3000)))...
   %          /sum(concat2.LTAS.norm0spectrum((concat2.LTAS.f>=50&concat2.LTAS.f<=10000)));
   %
   %       p = polyfit(concat2.LTAS.f(concat2.LTAS.f>50&concat2.LTAS.f<10000),...
   %          concat2.LTAS.dBspectrum(concat2.LTAS.f>50&concat2.LTAS.f<10000),1);
   %       concat2.LTAS.dBspect_slope=p;
   %       p = polyfit(concat2.LTAS.f(concat2.LTAS.f>Fo_stats.median&concat2.LTAS.f<10000),...
   %          concat2.LTAS.dBspectrum(concat2.LTAS.f>Fo_stats.median&concat2.LTAS.f<10000),1);
   %       concat2.LTAS.dBspect_slopeFo=p;
   %    end
   %    % semitone standard deviation
   %    concat2.STSD = 12*log2((concat2.Fo_stats.mean+concat2.Fo_stats.std/2)/(concat2.Fo_stats.mean-concat2.Fo_stats.std/2));
   %    concat2.STSDv = 12*log2((concat2.Fov_stats.mean+concat2.Fov_stats.std/2)/(concat2.Fov_stats.mean-concat2.Fov_stats.std/2));
   
   
   
   %         durOspch1: 20.7100
   %         durOspch2: 17.5600
   %         durOvoice: 17.1867
   %       durOvoicePr: 17.7100
   %        perOvoice1: 97.8740
   %        perOvoice2: 82.9873
   %     FoAllhz_stats: [1×1 struct]
   %     FoAllst_stats: [1×1 struct]
   %     FoPRThz_stats: [1×1 struct]
   %     FoPRTst_stats: [1×1 struct]
   %     FoSHRhz_stats: [1×1 struct]
   %     FoSHRst_stats: [1×1 struct]
   %     FoAUDhz_stats: [1×1 struct]
   %     FoAUDst_stats: [1×1 struct]
   %       SoAud_stats: [1×1 struct]
   %          dB_stats: [1×1 struct]
   %         dB_statsi: [1×1 struct]
   
   outdata=[];     % metrics
   outdataD=[];    % metrics names
   
   % Parameters
   outdata=[outdata,NaN,...
      u_trimYN,c_fryanaly,s_RMS_window,s_foS_step,...     %15
      s_foS_SilenceThrsh,s_foS_VoiceThrsh,s_foS_OctCost,s_foS_OctJumpCost,...
      s_foS_VoiceUnvoiceCost,s_pertS_MaxPeriodFact,s_pertS_MaxAmpFact,...
      NaN,NaN,c_LTAS_lower,c_LTAS_upper];
   outdataD=[outdataD;'*SETTINGS*';...
      'trimOption'; 'fry_analy_'; 'RMS_windw_'; 'Fo_WinStp_';...
      'perSilncTh'; 'perVoicTh_'; 'perOctCst_'; 'perOctJmp_';...
      'perVcUnVcC'; 'perMxPerF_'; 'perMxAmpF_'; ...
      '-spc-     '; '-spc-     '; 'LTASlwr_Bd'; 'LTASupprBd'];%slope accross utterance
   outdataD=[outdataD;'----------'];
   outdata=[outdata,NaN];
   
   outdata=[outdata,NaN,...      %
      concat2.durOspch1,...      % end voicing - beginning voicing
      concat2.durOspch2,...      % estimated length of all speaking parts
      concat2.durOvoice,...      % duration of voice only combining all 3 extraction types
      concat2.durOvoicePr,...    % duration of voice only from PRAAT raw
      concat2.perOvoice1,...     % voicing percentage compared to judged amount of speech (voicing time of the concatinated sections)
      concat2.perOvoice2,...     % voicing percentage compared to the beginning and end of the speech start (trimmed but with pauses)
      concat2.g_durationavg,...  % average length of read files
      concat2.g_durationtot];    % overall length of read files
   outdataD=[outdataD;'*GEN_DATA_';...
      'DurBeg2End';...      % end voicing - beginning voicing
      'DurSpchDtc';...      % estimated length of all speaking parts
      'DurVoicALL';...      % duration of voice only (fo segments from all 3 extraction methods)
      'DurVoicPRT';...      % duration of voice only (Praat only voicing)
      '%VoicCnSpc';...      % Voice % compared to the concatenated speech
      '%VoicAlSpc';...      % Voice % from the length of speech, from first instance to last instance
      'avgFileLen';...      % Voice % from the length of speech, from first instance to last instance
      'totFileLen'];
   outdataD=[outdataD;'----------'];
   outdata=[outdata,NaN];
   
   outdata=[outdata,NaN,...
      concat2.dB_stats.n,concat2.dB_stats.mean,concat2.dB_stats.std,concat2.dB_stats.var,...     %13
      concat2.dB_stats.mode,concat2.dB_stats.median,concat2.dB_stats.Q1,concat2.dB_stats.Q3,...
      concat2.dB_stats.IQR,concat2.dB_stats.skew,concat2.dB_stats.bowley_skew,concat2.dB_stats.kurt,...
      concat2.dB_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*dB_ALLraw';...
      'dB__N     '; 'dB__mean  '; 'dB__std   '; 'dB__var   ';...
      'dB__mode  '; 'dB__med   '; 'dB__Q1    '; 'dB__Q3    ';...
      'dB__IQR   '; 'dB__skew  '; 'dB__bowskw'; 'dB__kurt  ';...
      'dB__slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.dBv_stats.n,concat2.dBv_stats.mean,concat2.dBv_stats.std,concat2.dBv_stats.var,...     %13
      concat2.dBv_stats.mode,concat2.dBv_stats.median,concat2.dBv_stats.Q1,concat2.dBv_stats.Q3,...
      concat2.dBv_stats.IQR,concat2.dBv_stats.skew,concat2.dBv_stats.bowley_skew,concat2.dBv_stats.kurt,...
      concat2.dBv_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*dBvALLspc';...
      'dBv_N     '; 'dBv_mean  '; 'dBv_std   '; 'dBv_var   ';...
      'dBv_mode  '; 'dBv_med   '; 'dBv_Q1    '; 'dBv_Q3    ';...
      'dBv_IQR   '; 'dBv_skew  '; 'dBv_bowskw'; 'dBv_kurt  ';...
      'dBv_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.SoAud_stats.n,concat2.SoAud_stats.mean,concat2.SoAud_stats.std,concat2.SoAud_stats.var,...     %13
      concat2.SoAud_stats.mode,concat2.SoAud_stats.median,concat2.SoAud_stats.Q1,concat2.SoAud_stats.Q3,...
      concat2.SoAud_stats.IQR,concat2.SoAud_stats.skew,concat2.SoAud_stats.bowley_skew,concat2.SoAud_stats.kurt,...
      concat2.SoAud_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Sov_PchSt';...
      'Sov_N     '; 'Sov_mean  '; 'Sov_std   '; 'Sov_var   ';...
      'Sov_mode  '; 'Sov_med   '; 'Sov_Q1    '; 'Sov_Q3    ';...
      'Sov_IQR   '; 'Sov_skew  '; 'Sov_bowskw'; 'Sov_kurt  ';...
      'Sov_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoAllhz_stats.n,concat2.FoAllhz_stats.mean,concat2.FoAllhz_stats.std,concat2.FoAllhz_stats.var,...     %13
      concat2.FoAllhz_stats.mode,concat2.FoAllhz_stats.median,concat2.FoAllhz_stats.Q1,concat2.FoAllhz_stats.Q3,...
      concat2.FoAllhz_stats.IQR,concat2.FoAllhz_stats.skew,concat2.FoAllhz_stats.bowley_skew,concat2.FoAllhz_stats.kurt,...
      concat2.FoAllhz_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Fo_Hz_ALL';...
      'FoAl_N    '; 'FoAl_mean '; 'FoAl_std  '; 'FoAl_var  ';...
      'FoAl_mode '; 'FoAl_med  '; 'FoAl_Q1   '; 'FoAl_Q3   ';...
      'FoAl_IQR  '; 'FoAl_skew '; 'FoAlbowskw'; 'FoAl_kurt ';...
      'FoAl_slop '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoAllst_stats.n,concat2.FoAllst_stats.mean,concat2.FoAllst_stats.std,concat2.FoAllst_stats.var,...     %13
      concat2.FoAllst_stats.mode,concat2.FoAllst_stats.median,concat2.FoAllst_stats.Q1,concat2.FoAllst_stats.Q3,...
      concat2.FoAllst_stats.IQR,concat2.FoAllst_stats.skew,concat2.FoAllst_stats.bowley_skew,concat2.FoAllst_stats.kurt,...
      concat2.FoAllst_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Fo_ST_ALL';...
      'FoAl_N    '; 'FoAl_mean '; 'FoAl_std  '; 'FoAl_var  ';...
      'FoAl_mode '; 'FoAl_med  '; 'FoAl_Q1   '; 'FoAl_Q3   ';...
      'FoAl_IQR  '; 'FoAl_skew '; 'FoAlbowskw'; 'FoAl_kurt ';...
      'FoAl_slop '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoPRThz_stats.n,concat2.FoPRThz_stats.mean,concat2.FoPRThz_stats.std,concat2.FoPRThz_stats.var,...     %13
      concat2.FoPRThz_stats.mode,concat2.FoPRThz_stats.median,concat2.FoPRThz_stats.Q1,concat2.FoPRThz_stats.Q3,...
      concat2.FoPRThz_stats.IQR,concat2.FoPRThz_stats.skew,concat2.FoPRThz_stats.bowley_skew,concat2.FoPRThz_stats.kurt,...
      concat2.FoPRThz_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Fov_Hz___';...
      'Fov_N     '; 'Fov_mean  '; 'Fov_std   '; 'Fov_var   ';...
      'Fov_mode  '; 'Fov_med   '; 'Fov_Q1    '; 'Fov_Q3    ';...
      'Fov_IQR   '; 'Fov_skew  '; 'Fov_bowskw'; 'Fov_kurt  ';...
      'Fov_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoPRTst_stats.n,concat2.FoPRTst_stats.mean,concat2.FoPRTst_stats.std,concat2.FoPRTst_stats.var,...     %13
      concat2.FoPRTst_stats.mode,concat2.FoPRTst_stats.median,concat2.FoPRTst_stats.Q1,concat2.FoPRTst_stats.Q3,...
      concat2.FoPRTst_stats.IQR,concat2.FoPRTst_stats.skew,concat2.FoPRTst_stats.bowley_skew,concat2.FoPRTst_stats.kurt,...
      concat2.FoPRTst_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Fov_ST___';...
      'Fov_N     '; 'Fov_mean  '; 'Fov_std   '; 'Fov_var   ';...
      'Fov_mode  '; 'Fov_med   '; 'Fov_Q1    '; 'Fov_Q3    ';...
      'Fov_IQR   '; 'Fov_skew  '; 'Fov_bowskw'; 'Fov_kurt  ';...
      'Fov_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoAUDhz_stats.n,concat2.FoAUDhz_stats.mean,concat2.FoAUDhz_stats.std,concat2.FoAUDhz_stats.var,...     %13
      concat2.FoAUDhz_stats.mode,concat2.FoAUDhz_stats.median,concat2.FoAUDhz_stats.Q1,concat2.FoAUDhz_stats.Q3,...
      concat2.FoAUDhz_stats.IQR,concat2.FoAUDhz_stats.skew,concat2.FoAUDhz_stats.bowley_skew,concat2.FoAUDhz_stats.kurt,...
      concat2.FoAUDhz_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Pov_Hz___';...
      'Pov_N     '; 'Pov_mean  '; 'Pov_std   '; 'Pov_var   ';...
      'Pov_mode  '; 'Pov_med   '; 'Pov_Q1    '; 'Pov_Q3    ';...
      'Pov_IQR   '; 'Pov_skew  '; 'Pov_bowskw'; 'Pov_kurt  ';...
      'Pov_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoAUDst_stats.n,concat2.FoAUDst_stats.mean,concat2.FoAUDst_stats.std,concat2.FoAUDst_stats.var,...     %13
      concat2.FoAUDst_stats.mode,concat2.FoAUDst_stats.median,concat2.FoAUDst_stats.Q1,concat2.FoAUDst_stats.Q3,...
      concat2.FoAUDst_stats.IQR,concat2.FoAUDst_stats.skew,concat2.FoAUDst_stats.bowley_skew,concat2.FoAUDst_stats.kurt,...
      concat2.FoAUDst_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Pov_ST___';...
      'Pov_N     '; 'Pov_mean  '; 'Pov_std   '; 'Pov_var   ';...
      'Pov_mode  '; 'Pov_med   '; 'Pov_Q1    '; 'Pov_Q3    ';...
      'Pov_IQR   '; 'Pov_skew  '; 'Pov_bowskw'; 'Pov_kurt  ';...
      'Pov_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoSHRhz_stats.n,concat2.FoSHRhz_stats.mean,concat2.FoSHRhz_stats.std,concat2.FoSHRhz_stats.var,...     %13
      concat2.FoSHRhz_stats.mode,concat2.FoSHRhz_stats.median,concat2.FoSHRhz_stats.Q1,concat2.FoSHRhz_stats.Q3,...
      concat2.FoSHRhz_stats.IQR,concat2.FoSHRhz_stats.skew,concat2.FoSHRhz_stats.bowley_skew,concat2.FoSHRhz_stats.kurt,...
      concat2.FoSHRhz_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Fo_SHR_Hz';...
      'FoS_N     '; 'FoS_mean  '; 'FoS_std   '; 'FoS_var   ';...
      'FoS_mode  '; 'FoS_med   '; 'FoS_Q1    '; 'FoS_Q3    ';...
      'FoS_IQR   '; 'FoS_skew  '; 'FoS_bowskw'; 'FoS_kurt  ';...
      'FoS_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   outdata=[outdata,NaN,...
      concat2.FoSHRst_stats.n,concat2.FoSHRst_stats.mean,concat2.FoSHRst_stats.std,concat2.FoSHRst_stats.var,...     %13
      concat2.FoSHRst_stats.mode,concat2.FoSHRst_stats.median,concat2.FoSHRst_stats.Q1,concat2.FoSHRst_stats.Q3,...
      concat2.FoSHRst_stats.IQR,concat2.FoSHRst_stats.skew,concat2.FoSHRst_stats.bowley_skew,concat2.FoSHRst_stats.kurt,...
      concat2.FoSHRst_stats.slope(1)];
   outdata=[outdata,NaN];
   outdataD=[outdataD;'*Fo_SHR_ST';...
      'FoS_N     '; 'FoS_mean  '; 'FoS_std   '; 'FoS_var   ';...
      'FoS_mode  '; 'FoS_med   '; 'FoS_Q1    '; 'FoS_Q3    ';...
      'FoS_IQR   '; 'FoS_skew  '; 'FoS_bowskw'; 'FoS_kurt  ';...
      'FoS_slop  '];%slope accross utterance
   outdataD=[outdataD;'----------'];
   
   
   %    if u_concatYN==2
   %       outdata=[outdata,...
   %          concat2.ASLdata.activityFactor,...
   %          concat2.LTAS.totdBfft,LTAS.totdBrms,...
   %          concat2.LTAS.AlphaRatio,LTAS.dBspect_slope,LTAS.dBspect_slopeFo,...
   %          concat2.LTAS.dB1kTo3k,...
   %          concat2.LTAS.dB1kTo3kN,...
   %          concat2.c_HFCC.HFCC_SDS,...
   %          concat2.c_HFCC.DHFCC_SDS,...
   %          concat2.c_HFCC.DDHFCC_SDS];
   %       outdata=[outdata,NaN];
   %       outdataD=[outdataD;...
   %          'ActyFact';...      % activity factor activityFactor [0 1]
   %          'spdB_fft'; 'spdB_rms'; ...
   %          'AlphaRto'; 'dBspcSlp'; 'dBspcSpF';... % alpha ratio, spec slope from 50hz, spec slope from Fo
   %          'dB1kTo3k';...      % energy from 1kHz to 3kHz normed to that range
   %          'dB1kT3kN';...      % percent energy from 1kHz to 3kHz compared to 50hz-10khz
   %          'HFCCSDS ';...      % HFCC HFCC standard deviation sum measure
   %          'DHFCCSDS';...      % delta HFCC features, standard deviation sum measure
   %          'DDHFCSDS'];        % delta delta HFCC features, standard deviation sum measure
   %       outdataD=[outdataD;'--------'];
   %    end
   %
   %    if u_concatYN==2
   %       outdata=[outdata,concat2.jitter1.jitter,concat2.jitter1.jitter_abs,concat2.jitter1.jitter_rap,concat2.jitter1.jitter_ddp,...
   %          concat2.jitter1.shimmer,concat2.jitter1.shimmer_db,concat2.jitter1.shimmer_apq3,concat2.jitter1.shimmer_apq5,...
   %          concat2.jitter1.shimmer_apq11,concat2.jitter1.shimmer_dda,...
   %          concat2.jitter1.nhr,concat2.jitter1.hnr,...
   %          concat2.jitter1.rpde,concat2.jitter1.dfa,concat2.jitter1.ppe,...
   %          concat2.jitter1.length];        %14
   %       outdata=[outdata,NaN];
   %       outdataD=[outdataD;...
   %          'jitter  ';'jit_abs ';'jit_rap ';'jit_ddp ';...
   %          'shimmer ';'shim_db ';'shmapq3 ';'shmapq5 ';...
   %          'shmapq11';'shim_dda';...
   %          'nhr     ';'hnr     ';...
   %          'rpde    ';'dfa     ';'ppe     ';...
   %          'purt dur'];
   %       outdataD=[outdataD;'--------'];
   %    else
   %    end
   
   disp(' --end-- Concatenated analysis')
   
   
   %% save matlab file of results
   
   % save(['FullResults_' titlename '_mat.mat'],'cnt','filename','outdataD','outdataN')
   %    save([g_titlename '_concat1' '_mat.mat'],'concat1')
   save([g_titlename '_all' '_mat.mat'])
   
   if s_name == 0
      workSheet = 'data';
      xlswrite([s_files '\' g_titlename '_ConcatResults' '.xlsx'],[{'Filenames '} cellstr(outdataD)'],workSheet,'A1');
      xlswrite([s_files '\' g_titlename '_ConcatResults' '.xlsx'],g_filename,workSheet,'A3');
      xlswrite([s_files '\' g_titlename '_ConcatResults' '.xlsx'],outdata,workSheet,'B2');
   elseif s_name == 1
      xlswrite([s_files '\' g_outname '_ConcatResults' '.xlsx'],[{'Filenames'} cellstr(outdataD)'],workSheet,'A1');
      xlswrite([s_files '\' g_outname '_ConcatResults' '.xlsx'],g_filename,workSheet,'A3');
      xlswrite([s_files '\' g_outname '_ConcatResults' '.xlsx'],outdata,workSheet,'B2');
   end
   
   clear worksheet
   
   %% plot concatinated results
   if u_figson==1||u_figson==3
      %%
      fig1=6;fig2=7;
      figure(fig1),clf,subplot(4,1,1)
      rngtmp=[min(cncat1_i_t) max(cncat1_i_t)];
      tmp3=concat2.dBv_stats.slope(1)*rngtmp+concat2.dBv_stats.slope(2);
      tmp4=concat2.dB_stats.slope(1)*rngtmp+concat2.dB_stats.slope(2);
      plot(cncat1_i_t,cncat1_sig_dBi,'.r',rngtmp,tmp3,'k',rngtmp,tmp4,'k:',...
         rngtmp,[1 1]*concat2.dBv_stats.mean-concat2.dBv_stats.std,'b--',...
         rngtmp,[1 1]*concat2.dBv_stats.mean+concat2.dBv_stats.std,'b--')
      axis([0 max(cncat1_i_t) 0.95*min(cncat1_sig_dBi) 1.05*max(cncat1_sig_dBi) ])
      ylabel('dB'),
      tmp=strrep(cd,'_', ' ');tmp=strrep(tmp,'\','.');title(tmp,'FontSize',6)
      text(0.5*max(cncat1_i_t),0.95*max(cncat1_sig_dBi),['|slopeVoic=' num2str(concat2.dBv_stats.slope(1))],'FontSize',6)
      text(0.5*max(cncat1_i_t),0.80*max(cncat1_sig_dBi),[':slopeAll=' num2str(concat2.dB_stats.slope(1))],'FontSize',6)
      
      figure(fig1),subplot(4,1,2),
      rngtmp=[min(cncat1_t_t) max(cncat1_i_t)];
      tmp3=concat2.SoAud_stats.slope(1)*rngtmp+concat2.SoAud_stats.slope(2);
      plot(cncat1_t_t,cncat1_f0_So_aud,'.g',cncat1_i_t,cncat1_f0_So_audi,'.r',...
         rngtmp,tmp3,'m',rngtmp,[1 1]*concat2.SoAud_stats.mean,'b-',...
         rngtmp,[1 1]*(concat2.SoAud_stats.mean+concat2.SoAud_stats.std),'g--',...
         rngtmp,[1 1]*(concat2.SoAud_stats.mean-concat2.SoAud_stats.std),'b--')
      axis([0 max(cncat1_i_t) 0.96*concat2.SoAud_stats.min concat2.SoAud_stats.max])
      ylabel('P_s'),
      title(['slope=' num2str(round(10*concat2.SoAud_stats.slope(1))/10)],'Fontsize',6)
      
      figure(fig1),subplot(4,1,3),
      rngtmp=[min(cncat1_i_t) max(cncat1_i_t)];
      tmp3=concat2.FoAllhz_stats.slope(1)*rngtmp+concat2.FoAllhz_stats.slope(2);
      tmp3a=concat2.FoPRThz_stats.slope(1)*rngtmp+concat2.FoPRThz_stats.slope(2);
      tmp3b=concat2.FoAUDhz_stats.slope(1)*rngtmp+concat2.FoAUDhz_stats.slope(2);
      tmp3c=concat2.FoSHRhz_stats.slope(1)*rngtmp+concat2.FoSHRhz_stats.slope(2);
      plot(cncat1_i_t,cncat1_f0_fo_prtfi,'.k',cncat1_i_t,cncat1_f0_fo_audi,'.r',...
         cncat1_i_t,cncat1_f0_fo_shpi,'.b',rngtmp,tmp3,'m',rngtmp,tmp3a,'k:',...
         rngtmp,tmp3b,'r:',rngtmp,tmp3c,'b:',rngtmp,[1 1]*concat2.FoAllhz_stats.stmean,'b-',...
         rngtmp,[1 1]*concat2.FoAllhz_stats.stmeanmstd,'b--',...
         rngtmp,[1 1]*concat2.FoAllhz_stats.stmeanpstd,'b--',...
         rngtmp,[1 1]*concat2.FoAllhz_stats.stmeanp2std,'b:')
      axis([0 max(cncat1_i_t) 0.96*concat2.FoAllhz_stats.min concat2.FoAllhz_stats.max])
      ylabel('Hz'),
      title(['slope: allm=' num2str(round(10*concat2.FoAllhz_stats.slope(1))/10) ', prtk=' num2str(round(10*concat2.FoPRThz_stats.slope(1))/10) ...
         ', audr=' num2str(round(10*concat2.FoAUDhz_stats.slope(1))/10) ', shrb=' num2str(round(10*concat2.FoSHRhz_stats.slope(1))/10)],'Fontsize',6)
      
      figure(fig1),subplot(4,1,4),
      rngtmp=[min(cncat1_i_t) max(cncat1_i_t)];
      tmp3=concat2.FoAllst_stats.slope(1)*rngtmp+concat2.FoAllst_stats.slope(2);
      tmp3a=concat2.FoPRTst_stats.slope(1)*rngtmp+concat2.FoPRTst_stats.slope(2);
      tmp3b=concat2.FoAUDst_stats.slope(1)*rngtmp+concat2.FoAUDst_stats.slope(2);
      tmp3c=concat2.FoSHRst_stats.slope(1)*rngtmp+concat2.FoSHRst_stats.slope(2);
      plot(cncat1_i_t,12*log2(double(cncat1_f0_fo_prtfi)./220),'.k',...
         cncat1_i_t,12*log2(double(cncat1_f0_fo_audi)./220),'.r',...
         cncat1_i_t,12*log2(double(cncat1_f0_fo_shpi)./220),'.b',...
         rngtmp,tmp3,'m',rngtmp,tmp3a,'k:',rngtmp,tmp3b,'r:',rngtmp,tmp3c,'b:',...
         rngtmp,[1 1]*12*log2(concat2.FoAllhz_stats.stmean/220),'b-',...
         rngtmp,[1 1]*12*log2(concat2.FoAllhz_stats.stmeanmstd/220),'b--',...
         rngtmp,[1 1]*12*log2(concat2.FoAllhz_stats.stmeanpstd/220),'b--',...
         rngtmp,[1 1]*12*log2(concat2.FoAllhz_stats.stmeanp2std/220),'b:')
      axis([0 max(cncat1_i_t) 1.01*concat2.FoAllst_stats.min concat2.FoAllst_stats.max])
      ylabel('ST'),xlabel('Time (min)')
      title(['slope: allm=' num2str(round(10*concat2.FoAllst_stats.slope(1))/10) ', prtk=' num2str(round(10*concat2.FoPRTst_stats.slope(1))/10) ...
         ', audr=' num2str(round(10*concat2.FoAUDst_stats.slope(1))/10) ', shrb=' num2str(round(10*concat2.FoSHRst_stats.slope(1))/10)],'Fontsize',6)
      
      figure(fig1)
      saveas(gcf,['1' g_titlename '_Overlap' '.fig'])
      saveas(gcf,['1' g_titlename '_Overlap' '.tif'],'tiff')
      
      figure(fig2)
      % Fo (hz)
      subplot(2,2,1)
      rngtmp=[(concat2.FoAllhz_stats.stmeanm2std) (concat2.FoAllhz_stats.stmeanp2std)];
      rngtmp=[floor(rngtmp(1)) ceil(rngtmp(2))];
      rngtmp=[rngtmp(1):floor((rngtmp(2)-rngtmp(1))/20):rngtmp(2)];
      if rngtmp(2)-rngtmp(1)<20
         tmp=rngtmp(2)-rngtmp(1); tmp=ceil((20-tmp)/2);
         rngtmp(1)=rngtmp(1)-tmp; rngtmp(2)=rngtmp(2)+tmp;
      end      
      [tmpy,tmpx]=hist([cncat1_f0_fo_prtfi cncat1_f0_fo_audi cncat1_f0_fo_shpi],rngtmp);
      tmpy=100*tmpy/sum(tmpy);
      [tmpya,tmpxa]=hist([cncat1_f0_fo_prtfi],rngtmp);tmpya=100*tmpya/sum(tmpya);
      [tmpyb,tmpxb]=hist([cncat1_f0_fo_audi],rngtmp);tmpyb=100*tmpyb/sum(tmpyb);
      [tmpyc,tmpxc]=hist([cncat1_f0_fo_shpi],rngtmp);tmpyc=100*tmpyc/sum(tmpyc);
      tmpmax=1.1*max([tmpy tmpya tmpyb tmpyc]);
      plot(tmpx,tmpy,'.-',tmpxa,tmpya,'k.-',tmpxb,tmpyb,'b.--',tmpxc,tmpyc,'r.:',...
         [1 1]*concat2.FoAllhz_stats.stmean,[1 1.1*tmpmax],'g-',...
         [1 1]*concat2.FoAllhz_stats.stmeanmstd,[1 tmpmax],'g--',...
         [1 1]*concat2.FoAllhz_stats.stmeanpstd,[1 tmpmax],'g--',...
         [1 1]*concat2.FoAllhz_stats.median,[1 1.1*tmpmax],'m-',...
         [1 1]*concat2.FoAllhz_stats.Q1,[1 tmpmax],'m--',...
         [1 1]*concat2.FoAllhz_stats.Q3,[1 tmpmax],'m--')
      text(double(0.95*concat2.FoAllhz_stats.median),1.05*tmpmax,'med','Fontsize',7)
      text(double(0.95*concat2.FoAllhz_stats.Q1),tmpmax,'Q1','Fontsize',6)
      text(double(0.95*concat2.FoAllhz_stats.Q3),tmpmax,'Q3','Fontsize',6)
      text(double(0.95*concat2.FoAllhz_stats.stmean),0.96*tmpmax,'st\mu','Fontsize',6)
      title(['med=' num2str(round(10*concat2.FoAllhz_stats.median)/10) ...
         'hz; IQR=' num2str(round(10*concat2.FoAllhz_stats.IQR)/10) 'hz'],'Fontsize',7)
      xlabel('f_o (hz)')
      
      % Fo (ST)
      subplot(2,2,2)
      rngtmp=[(concat2.FoAllhz_stats.stmeanm2std) (concat2.FoAllhz_stats.stmeanp2std)];
      rngtmp=12*log2(rngtmp/220);rngtmp=[floor(rngtmp(1)) ceil(rngtmp(2))];
      rngtmp=[rngtmp(1):max([floor((rngtmp(2)-rngtmp(1))/20) 1]):rngtmp(2)];
      [tmpy,tmpx]=hist(12*log2([cncat1_f0_fo_prtfi cncat1_f0_fo_audi cncat1_f0_fo_shpi]/220),rngtmp);
      tmpy=100*tmpy/sum(tmpy);
      [tmpya,tmpxa]=hist(12*log2(cncat1_f0_fo_prtfi/220),rngtmp);tmpya=100*tmpya/sum(tmpya);
      [tmpyb,tmpxb]=hist(12*log2(cncat1_f0_fo_audi/220),rngtmp);tmpyb=100*tmpyb/sum(tmpyb);
      [tmpyc,tmpxc]=hist(12*log2(cncat1_f0_fo_shpi/220),rngtmp);tmpyc=100*tmpyc/sum(tmpyc);
      tmpmax=1.1*max([tmpy tmpya tmpyb tmpyc]);
      plot(tmpx,tmpy,'.-',tmpxa,tmpya,'k.-',tmpxb,tmpyb,'b.--',tmpxc,tmpyc,'r.:',...
         [1 1]*12*log2(concat2.FoAllhz_stats.stmean/220),[1 1.1*tmpmax],'m--',...
         [1 1]*12*log2(concat2.FoAllhz_stats.stmeanmstd/220),[1 tmpmax],'m:',...
         [1 1]*12*log2(concat2.FoAllhz_stats.stmeanpstd/220),[1 tmpmax],'m:',...
         [1 1]*12*log2(concat2.FoAllhz_stats.median/220),[1 1.1*tmpmax],'g',...
         [1 1]*12*log2(concat2.FoAllhz_stats.Q1/220),[1 tmpmax],'g--',...
         [1 1]*12*log2(concat2.FoAllhz_stats.Q3/220),[1 tmpmax],'g--')
      text(double(12*log2(.93*concat2.FoAllhz_stats.median/220)),0.95*tmpmax,'med','Fontsize',6)
      text(double(12*log2(.96*concat2.FoAllhz_stats.stmean/220)),1.05*tmpmax,'\mu','Fontsize',8)
      text(double(12*log2(.94*concat2.FoAllhz_stats.stmeanmstd/220)),tmpmax,'\mu-\sigma','Fontsize',7)
      text(double(12*log2(.94*concat2.FoAllhz_stats.stmeanpstd/220)),tmpmax,'\mu+\sigma','Fontsize',7)
      title(['ST \mu=' num2str(round(10*concat2.FoAllhz_stats.stmean)/10) ...
         'hz; ST \mu\pm\sigma=' num2str(round(10*concat2.FoAllhz_stats.stmeanpstd)/10) ' & '...
         num2str(round(10*concat2.FoAllhz_stats.stmeanmstd)/10) 'hz'],'Fontsize',6)
      xlabel('f_o (st)')
      
      % P S
      subplot(2,2,3)
      %       rngtmp=[concat2.SoAud_stats.mean-3*concat2.SoAud_stats.std concat2.SoAud_stats.mean+3*concat2.SoAud_stats.std];
      rngtmp=[0 75]; rngtmp=[rngtmp(1):floor((rngtmp(2)-rngtmp(1))/25):rngtmp(2)];
      [tmpy,tmpx]=hist([cncat1_f0_So_audi],rngtmp);
      tmpy=100*tmpy/sum(tmpy); tmpmax=max(tmpy);
      plot(tmpx,tmpy,'.-',...
         [1 1]*concat2.SoAud_stats.median,[1 1.1*tmpmax],'m-',...
         [1 1]*concat2.SoAud_stats.Q1,[1 tmpmax],'m--',...
         [1 1]*concat2.SoAud_stats.Q3,[1 tmpmax],'m--',...
         [1 1]*concat2.SoAud_stats.mean,[1 1.1*tmpmax],'b--',...
         [1 1]*(concat2.SoAud_stats.mean+concat2.SoAud_stats.std),[1 tmpmax],'b:',...
         [1 1]*(concat2.SoAud_stats.mean-concat2.SoAud_stats.std),[1 tmpmax],'b:')
      text(double(0.95*concat2.SoAud_stats.median),1.1*tmpmax,'med','Fontsize',6)
      text(double(0.95*concat2.SoAud_stats.Q1),1.05*tmpmax,'Q1','Fontsize',6)
      text(double(0.95*concat2.SoAud_stats.Q3),1.05*tmpmax,'Q3','Fontsize',6)
      text(double(0.95*concat2.SoAud_stats.mean),1*tmpmax,'\mu','Fontsize',6)
      text(double(0.95*concat2.SoAud_stats.mean+concat2.SoAud_stats.std),.9*tmpmax,'\mu+\sigma','Fontsize',6)
      text(double(0.95*concat2.SoAud_stats.mean-concat2.SoAud_stats.std),.9*tmpmax,'\mu-\sigma','Fontsize',6)
      title(['med=' num2str(round(10*concat2.SoAud_stats.median)/10) ...
         '; IQR=' num2str(round(10*concat2.SoAud_stats.IQR)/10) ...
         '; \mu=' num2str(round(10*concat2.SoAud_stats.mean)/10) ...
         '; \sigma=' num2str(round(10*concat2.SoAud_stats.std)/10)],'Fontsize',7)
      xlabel('P_s')
      
      % dBv
      subplot(2,2,4)
      rngtmp=[concat2.dBv_stats.min concat2.dBv_stats.max];
      rngtmp2=[concat2.dB_stats.min concat2.dB_stats.max];
      rngtmp2=floor(rngtmp2(1)):max([floor((rngtmp(2)-rngtmp(1))/25) 1]):ceil(rngtmp2(2));
      rngtmp=floor(rngtmp(1)):max([floor((rngtmp(2)-rngtmp(1))/25) 1]):ceil(rngtmp(2));
      [tmpy,tmpx]=hist([cncat1_sig_dBi],rngtmp);
      [tmpya,tmpxa]=hist([cncat1_sig_dB],rngtmp2);
      tmpy=100*tmpy/sum(tmpy); tmpya=max(tmpy)*tmpya/max(tmpya); tmpmax=max([tmpy]);
      
      try dB_clust=findclusters(cncat1_sig_dB'); catch, dB_clust=findclusters(cncat1_sig_dB(cncat1_sig_dB>1)');end
      sig_dB_noise_dist = pdf(dB_clust.low.gcluster,rngtmp2); %sig_dB_noise_dist=sig_dB_noise_dist/max(sig_dB_noise_dist);
      sig_dB_speech_dist = pdf(dB_clust.high.gcluster,rngtmp2); sig_dB_speech_dist=sig_dB_speech_dist/max(sig_dB_speech_dist);
      [~, tmpi] = min(abs(sig_dB_noise_dist./sig_dB_speech_dist-1)); %find intersection of the two pdf
      tmp=dB_clust.high.stats.mean-2.7*dB_clust.high.stats.std;
      dB_thresh=mean([rngtmp2(tmpi) tmp]); dB_thresh=round(dB_thresh*10)/10;
      
      plot(tmpxa,tmpya,'r.-',tmpx,tmpy,'b.-',[1 1]*dB_thresh,[0 1.2*tmpmax],...
         rngtmp2,3*tmpmax*sig_dB_noise_dist,'m--',rngtmp2,tmpmax*sig_dB_speech_dist,'b:',...
         [1 1]*concat2.dBv_stats.median,[1 1.2*tmpmax],'m-',...
         [1 1]*concat2.dBv_stats.Q1,[1 tmpmax],'m--',...
         [1 1]*concat2.dBv_stats.Q3,[1 tmpmax],'m--',...
         [1 1]*concat2.dBv_stats.mean,[1 1.2*tmpmax],'b--',...
         [1 1]*(concat2.dBv_stats.mean+concat2.dBv_stats.std),[1 tmpmax],'b:',...
         [1 1]*(concat2.dBv_stats.mean-concat2.dBv_stats.std),[1 tmpmax],'b:')
      axis([rngtmp2(1) rngtmp2(end) 0 1.3*tmpmax])
      text(double(0.98*concat2.dBv_stats.median),1.1*tmpmax,'med','Fontsize',6)
      text(double(0.98*concat2.dBv_stats.Q1),1.05*tmpmax,'Q1','Fontsize',6)
      text(double(0.98*concat2.dBv_stats.Q3),1.05*tmpmax,'Q3','Fontsize',6)
      text(double(0.98*concat2.dBv_stats.mean),1*tmpmax,'\mu','Fontsize',6)
      text(double(0.98*concat2.dBv_stats.mean+concat2.dBv_stats.std),.9*tmpmax,'\mu+\sigma','Fontsize',6)
      text(double(0.98*concat2.dBv_stats.mean-concat2.dBv_stats.std),.9*tmpmax,'\mu-\sigma','Fontsize',6)
      title(['IQR=' num2str(round(10*concat2.dBv_stats.IQR)/10) ...
         'dB; \sigma=' num2str(round(10*concat2.dBv_stats.std)/10) 'dB'],'Fontsize',7)
      xlabel('dB')
      
      figure(fig2)
      saveas(gcf,['1' g_titlename '_Histogram' '.fig'])
      saveas(gcf,['1' g_titlename '_Histogram' '.tif'],'tiff')
      %         figure(4)
      %         saveas(gcf,[ fname '_' titlename '_hist' '.fig'])
      %         saveas(gcf,[ fname '_' titlename '_hist' '.tif'],'tiff')
      
      clear gcluster1 gcluster2 idx obj p tmp3a tmp3b tmp3c tmpi tmpmax
      clear tmpt tmpt2 tmpxa tmpxb tmpxc tmpya tmpyb tmpyc fig2 fig1
      clear tmpstr tmp tmp1 tmp2 tmp3 tmp4 tmp1x tmp1y tmp2 tmp2x tmp2y tmpx tmpy tmpstats tmpstatus
      clear tmp_cluster_lower tmp_cluster_lowerg tmp_cluster_lower_dist tmp_cluster_lower_stats
      clear tmp_cluster_upper tmp_cluster_upperg tmp_cluster_upper_dist tmp_cluster_upper_stats
   end
   
end

%% failed file list

g_fid_out = fopen([g_titlename '_trackcomplete.txt']);
badfiles = load([g_titlename '_trackcomplete.txt']);
fclose(g_fid_out); %close file

if ~isempty(badfiles)
   tmpstr=[g_titlename '_SKIPPED_FILES'];
   dos(['mkdir ' tmpstr],'-echo');
   badfiles=-badfiles((badfiles<0));
   fprintf('\n %d files were skipped:\n',length(badfiles))
   for n=1:length(badfiles)
      fprintf('%s\n',(g_filename{badfiles(n)}))
      tmpstatus=copyfile(g_filename{badfiles(n)}, tmpstr);      
   end
end
%% move figure files into new folder
tmpstr=[g_titlename '_PLOTS'];
dos(['mkdir ' tmpstr],'-echo');

tmpstatus=movefile('_fig_*.fig', tmpstr);
if tmpstatus==0
   disp(' problem moving figures (.fig files) ')
end
tmpstatus=movefile('_tif_*.tif', tmpstr);
if tmpstatus==0
   disp(' problem moving figures (.tif files) ')
end


%% move results files into new folder
tmpstr=[g_titlename '_RESULTS'];
dos(['mkdir ' tmpstr],'-echo');

tmpstatus=movefile([g_titlename '_FullResults.txt'], tmpstr);
if tmpstatus==0
   disp(' problem moving results.txt file ')
end
tmpstatus=movefile([g_titlename '*_mat.mat'], tmpstr);
if tmpstatus==0
   disp(' problem moving results.mat file(s) ')
end
tmpstatus=movefile([g_titlename '*.xlsx'], tmpstr);
if tmpstatus==0
   disp(' problem moving results.xlsx file(s) ')
end
tmpstatus=movefile([g_titlename '_contour.csv'], tmpstr);
if tmpstatus==0
   disp(' problem moving countour.csv file(s) ')
end
tmpstatus=movefile(['1' g_titlename '_Histogram.fig'], tmpstr);
tmpstatus=movefile(['1' g_titlename '_Histogram.tif'], tmpstr)+tmpstatus;
tmpstatus=movefile(['1' g_titlename '_Overlap.fig'], tmpstr)+tmpstatus;
tmpstatus=movefile(['1' g_titlename '_Overlap.tif'], tmpstr)+tmpstatus;

if tmpstatus==0
   disp(' problem moving result fig file(s) ')
end

tmpstatus=copyfile([mfilename '.m'], tmpstr);
if tmpstatus==0
   disp(' problem copying analysis file ')
end


%% ending

disp(' '),disp('------------------------------------------- ')
disp(['  Finished ' num2str(g_cnt) ' files' ])
disp('------------------------------------------- '),disp(' ')
toc
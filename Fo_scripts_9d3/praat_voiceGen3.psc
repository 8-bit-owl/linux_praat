form Voice_analysis
   text infile ''
   text outfile ''
   real fo_lower 0
   real fo_upper 0
   real fo_NCand 0 
   text fo_Accuracy '' 
   real fo_SilenceThrsh 0 
   real fo_VoiceThrsh 0 
   real fo_OctCost 0 
   real fo_OctJumpCost 0 
   real fo_VoiceUnvoiceCost 0 
   real fo_step 0
   text analysistype ''
   real maxPeriodFact 0 
   real maxAmpFact 0 
endform
Read from file... 'infile$'
infile$ = selected$("Sound")

To Pitch ('analysistype$')... fo_step fo_lower fo_NCand 'fo_Accuracy$' fo_SilenceThrsh fo_VoiceThrsh fo_OctCost fo_OctJumpCost fo_VoiceUnvoiceCost fo_upper
select Sound 'infile$'
To PointProcess (periodic, 'analysistype$')... fo_lower fo_upper

select PointProcess 'infile$'
plus Pitch 'infile$'
plus Sound 'infile$'
voiceReport$ = Voice report... 0 0 fo_lower fo_upper maxPeriodFact maxAmpFact fo_SilenceThrsh fo_VoiceThrsh

jitter = extractNumber (voiceReport$, "Jitter (local): ")
jitter_abs = extractNumber (voiceReport$, "Jitter (local, absolute): ")
jitter_rap = extractNumber (voiceReport$, "Jitter (rap): ")
jitter_ppq5 = extractNumber (voiceReport$, "Jitter (ppq5): ")
jitter_ddp = extractNumber (voiceReport$, "Jitter (ddp): ")

shimmer = extractNumber (voiceReport$, "Shimmer (local): ")
shimmer_db = extractNumber (voiceReport$, "Shimmer (local, dB): ")
shimmer_apq3 = extractNumber (voiceReport$, "Shimmer (apq3): ")
shimmer_apq5 = extractNumber (voiceReport$, "Shimmer (apq5): ")
shimmer_apq11 = extractNumber (voiceReport$, "Shimmer (apq11): ")
shimmer_dda = extractNumber (voiceReport$, "Shimmer (dda): ")

nhr = extractNumber (voiceReport$, "Mean noise-to-harmonics ratio: ")
hnr = extractNumber (voiceReport$, "Mean harmonics-to-noise ratio: ")

filedelete 'outfile$'
fileappend 'outfile$' 'jitter:8''newline$'
fileappend 'outfile$' 'jitter_abs:8''newline$'
fileappend 'outfile$' 'jitter_rap:8''newline$'
fileappend 'outfile$' 'jitter_ppq5:8''newline$'
fileappend 'outfile$' 'jitter_ddp:8''newline$'
fileappend 'outfile$' 'shimmer:8''newline$'
fileappend 'outfile$' 'shimmer_db:8''newline$'
fileappend 'outfile$' 'shimmer_apq3:8''newline$'
fileappend 'outfile$' 'shimmer_apq5:8''newline$'
fileappend 'outfile$' 'shimmer_apq11:8''newline$'
fileappend 'outfile$' 'shimmer_dda:8''newline$'
fileappend 'outfile$' 'nhr:8''newline$'
fileappend 'outfile$' 'hnr:8''newline$'

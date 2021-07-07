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
endform
Read from file... 'infile$'
To Pitch ('analysistype$')... fo_step fo_lower fo_NCand 'fo_Accuracy$' fo_SilenceThrsh fo_VoiceThrsh fo_OctCost fo_OctJumpCost fo_VoiceUnvoiceCost fo_upper
Down to PitchTier
Write to headerless spreadsheet file... 'outfile$'
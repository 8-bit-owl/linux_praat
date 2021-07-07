form Voice_analysis
   text name ''
endform
Read from file... 'name$'
To Pitch (pitchAnalysis)... fo_step fo_lower pitchNCand pitchAccuracy pitchSilenceThrsh pitchVoiceThrsh pitchOctCost pitchOctJumpCost pitchVoiceUnvoiceCost fo_upper
Down to PitchTier
Write to headerless spreadsheet file... pitch.txt
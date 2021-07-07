form Voice_analysis
   text name ''
endform
Read from file... 'name$'
To Pitch (analysistype)... 0.01 fo_lower 15 yes 0.03 VoiceThresh 0.0025 0.35 0.20 fo_upper
Down to PitchTier
Write to headerless spreadsheet file... pitch.txt
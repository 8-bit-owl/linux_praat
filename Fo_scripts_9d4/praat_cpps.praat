form Voice_analysis
   text name ''
endform
Read from file... 'name$'
name$ = selected$("Sound")

durationAll = Get total duration
minimumSPL = Get minimum... 0 0 None
maximumSPL = Get maximum... 0 0 None

# Narrow-band spectrogram and LTAS

To Spectrogram... 0.03 4000 0.002 20 Gaussian
select Sound 'name$'
To Ltas... 1
minimumSpectrum = Get minimum... 0 4000 None
maximumSpectrum = Get maximum... 0 4000 None

# Power-cepstrogram, Cepstral peak prominence and Smoothed cepstral peak prominence

select Sound 'name$'
To PowerCepstrogram... 60 0.002 5000 50
cpps = Get CPPS... no 0.01 0.001 60 330 0.05 Parabolic 0.001 0 Straight Robust
To PowerCepstrum (slice)... 0.1
maximumCepstrum = Get peak... 60 330 None


# Calculation of the AVQI

#avqi = cpps

filedelete avqi2.txt 
fileappend avqi2.txt 'cpps:8''newline$'






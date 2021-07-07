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

# Slope of the long-term average spectrum

select Sound 'name$'
To Ltas... 1
slope = Get slope... 0 1000 1000 10000 energy

# Tilt of trendline through the long-term average spectrum

select Ltas 'name$'
Compute trend line... 1 10000
tilt = Get slope... 0 1000 1000 10000 energy

# Amplitude perturbation measures

select Sound 'name$'
To PointProcess (periodic, cc)... 50 400
Rename... avqi1
select Sound 'name$'
plus PointProcess avqi1
percentShimmer = Get shimmer (local)... 0 0 0.0001 0.02 1.3 1.6
shim = percentShimmer*100
shdb = Get shimmer (local_dB)... 0 0 0.0001 0.02 1.3 1.6

# Harmonic-to-noise ratio

select Sound 'name$'
To Pitch (cc)... 0 75 15 no 0.03 0.45 0.01 0.35 0.14 600
select Sound 'name$'
plus Pitch 'name$'
To PointProcess (cc)
Rename... avqi2
select Sound 'name$'
plus Pitch 'name$'
plus PointProcess avqi2
voiceReport$ = Voice report... 0 0 75 600 1.3 1.6 0.03 0.45
hnr = extractNumber (voiceReport$, "Mean harmonics-to-noise ratio: ")

# Calculation of the AVQI

#avqi = ((3.295-(0.111*cpps)-(0.073*hnr)-(0.213*shim)+(2.789*shdb)-#(0.032*slope)+(0.077*tilt))*2.208)+1.797

filedelete avqi2.txt 
fileappend avqi2.txt 'cpps:8''newline$'
fileappend avqi2.txt 'hnr:8''newline$'
fileappend avqi2.txt 'shim:8''newline$'
fileappend avqi2.txt 'shdb:8''newline$'
fileappend avqi2.txt 'slope:8''newline$'
fileappend avqi2.txt 'tilt:8''newline$'
#fileappend avqi2.txt 'avqi:8''newline$'







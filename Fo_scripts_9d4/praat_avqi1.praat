# --------------------------------------------------------------------------------------------
# PART 0:
# HIGH-PASS FILTERING OF THE SOUND FILES. 
# --------------------------------------------------------------------------------------------

Read from file... cs.wav
Read from file... sv.wav
select Sound cs
Filter (stop Hann band)... 0 34 0.1
Rename... cs
select Sound sv
Filter (stop Hann band)... 0 34 0.1
Rename... sv

# --------------------------------------------------------------------------------------------
# PART 1:
# DETECTION, EXTRACTION AND CONCATENATION OF
# THE VOICED SEGMENTS IN THE RECORDING
# OF CONTINUOUS SPEECH.
# --------------------------------------------------------------------------------------------

select Sound cs
Copy... original
samplingRate = Get sampling frequency
intermediateSamples = Get sampling period
Create Sound... onlyVoice 0 0.001 'samplingRate' 0 
select Sound original
To TextGrid (silences)... 50 0.003 -25 0.1 0.1 silence sounding
select Sound original
plus TextGrid original
Extract intervals where... 1 no "does not contain" silence
Concatenate
select Sound chain
Rename... onlyLoud
globalPower = Get power in air
select TextGrid original
Remove

select Sound onlyLoud
signalEnd = Get end time
windowBorderLeft = Get start time
windowWidth = 0.03
windowBorderRight = windowBorderLeft + windowWidth
globalPower = Get power in air
voicelessThreshold = globalPower*(30/100)

select Sound onlyLoud
extremeRight = signalEnd - windowWidth
while windowBorderRight < extremeRight
	Extract part... 'windowBorderLeft' 'windowBorderRight' Rectangular 1.0 no
	select Sound onlyLoud_part
	partialPower = Get power in air
	if partialPower > voicelessThreshold
		call checkZeros 0
		if (zeroCrossingRate <> undefined) and (zeroCrossingRate < 3000)
			select Sound onlyVoice
			plus Sound onlyLoud_part
			Concatenate
			Rename... onlyVoiceNew
			select Sound onlyVoice
			Remove
			select Sound onlyVoiceNew
			Rename... onlyVoice
		endif
	endif
	select Sound onlyLoud_part
	Remove
	windowBorderLeft = windowBorderLeft + 0.03
	windowBorderRight = windowBorderLeft + 0.03
	select Sound onlyLoud
endwhile
select Sound onlyVoice

procedure checkZeros zeroCrossingRate

	start = 0.0025
	startZero = Get nearest zero crossing... 'start'
	findStart = startZero
	findStartZeroPlusOne = startZero + intermediateSamples
	startZeroPlusOne = Get nearest zero crossing... 'findStartZeroPlusOne'
	zeroCrossings = 0
	strips = 0

	while (findStart < 0.0275) and (findStart <> undefined)
		while startZeroPlusOne = findStart
			findStartZeroPlusOne = findStartZeroPlusOne + intermediateSamples
			startZeroPlusOne = Get nearest zero crossing... 'findStartZeroPlusOne'
		endwhile
		afstand = startZeroPlusOne - startZero
		strips = strips +1
		zeroCrossings = zeroCrossings +1
		findStart = startZeroPlusOne
	endwhile
	zeroCrossingRate = zeroCrossings/afstand
endproc

# --------------------------------------------------------------------------------------------
# PART 2:
# DETERMINATION OF THE SIX ACOUSTIC MEASURES
# AND CALCULATION OF THE ACOUSTIC VOICE QUALITY INDEX.
# --------------------------------------------------------------------------------------------

select Sound sv
durationVowel = Get total duration
durationStart=durationVowel-3
if durationVowel>3
Extract part... durationStart durationVowel rectangular 1 no
Rename... sv2
elsif durationVowel<=3
Copy... sv2
endif

# test output
# select Sound sv2
# Save as WAV file: "svtest.wav"
# select Sound onlyVoice
# Save as WAV file: "OVtest.wav"

select Sound onlyVoice
durationOnlyVoice = Get total duration
plus Sound sv2
Concatenate
Rename... avqi
durationAll = Get total duration
minimumSPL = Get minimum... 0 0 None
maximumSPL = Get maximum... 0 0 None

# select Sound avqi
# Save as WAV file: "AVQItest.wav"

# Narrow-band spectrogram and LTAS

To Spectrogram... 0.03 4000 0.002 20 Gaussian
select Sound avqi
To Ltas... 1
minimumSpectrum = Get minimum... 0 4000 None
maximumSpectrum = Get maximum... 0 4000 None

# Power-cepstrogram, Cepstral peak prominence and Smoothed cepstral peak prominence

select Sound avqi
To PowerCepstrogram... 60 0.002 5000 50
cpps = Get CPPS... no 0.01 0.001 60 330 0.05 Parabolic 0.001 0 Straight Robust
To PowerCepstrum (slice)... 0.1
maximumCepstrum = Get peak... 60 330 None

# Slope of the long-term average spectrum

select Sound avqi
To Ltas... 1
slope = Get slope... 0 1000 1000 10000 energy

# Tilt of trendline through the long-term average spectrum

select Ltas avqi
Compute trend line... 1 10000
tilt = Get slope... 0 1000 1000 10000 energy

# Amplitude perturbation measures

select Sound avqi
To PointProcess (periodic, cc)... 50 400
Rename... avqi1
select Sound avqi
plus PointProcess avqi1
percentShimmer = Get shimmer (local)... 0 0 0.0001 0.02 1.3 1.6
shim = percentShimmer*100
shdb = Get shimmer (local_dB)... 0 0 0.0001 0.02 1.3 1.6

# Harmonic-to-noise ratio

select Sound avqi
To Pitch (cc)... 0 75 15 no 0.03 0.45 0.01 0.35 0.14 600
select Sound avqi
plus Pitch avqi
To PointProcess (cc)
Rename... avqi2
select Sound avqi
plus Pitch avqi
plus PointProcess avqi2
voiceReport$ = Voice report... 0 0 75 600 1.3 1.6 0.03 0.45
hnr = extractNumber (voiceReport$, "Mean harmonics-to-noise ratio: ")

# Calculation of the AVQI

avqi = ((3.295-(0.111*cpps)-(0.073*hnr)-(0.213*shim)+(2.789*shdb)-(0.032*slope)+(0.077*tilt))*2.208)+1.797

filedelete avqi1.txt 
fileappend avqi1.txt 'cpps:8''newline$'
fileappend avqi1.txt 'hnr:8''newline$'
fileappend avqi1.txt 'shim:8''newline$'
fileappend avqi1.txt 'shdb:8''newline$'
fileappend avqi1.txt 'slope:8''newline$'
fileappend avqi1.txt 'tilt:8''newline$'
fileappend avqi1.txt 'avqi:8''newline$'







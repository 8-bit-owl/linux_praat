form Voice_analysis
   text name ''
endform
Read from file... 'name$'
name$ = selected$("Sound")

select Sound 'name$'
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

select Sound onlyVoice
Save as WAV file: "pconcat.wav"
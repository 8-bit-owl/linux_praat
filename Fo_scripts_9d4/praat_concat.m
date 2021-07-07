function pconcat = praat_concat(sig,Fs)

audiowrite('toconcat.wav',sig,Fs)
[s w] = system('Praat.exe --run praat_concat.praat toconcat.wav');

delete('toconcat.wav')
[pconcat fs] = audioread('pconcat.wav');
delete('pconcat.wav')


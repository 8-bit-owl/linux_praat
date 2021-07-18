function pconcat = praat_concat(sig,Fs)

audiowrite('toconcat.wav',sig,Fs)
if isunix == 1; [s w] = system('Praat.exe --run praat_concat.praat toconcat.wav');
elseif ispc == 1; [s w] = system('Praat.exe --run praat_concat.praat toconcat.wav');
end
delete('toconcat.wav')
[pconcat fs] = audioread('pconcat.wav');
delete('pconcat.wav')


function outsig = avqi3(insig,Fs)

audiowrite('avqi2.wav',insig,Fs)

[s w] = system('Praat.exe --run praat_avqi2.praat avqi2.wav');
datatemp = importdata('avqi2.txt');
outsig.cpps = datatemp(1);
outsig.hnr = datatemp(2);
outsig.shim = datatemp(3);
outsig.shdb = datatemp(4);
outsig.slope = datatemp(5);
outsig.tilt = datatemp(6);
outsig.avqi = NaN;
delete('avqi2.txt')
delete('avqi2.wav')

function outavqi2 = avqi2(concatsig,svsig,fs)

if length(svsig)/fs > 3
    svsig = svsig(end-3*fs:end);
end

avqiwav2 = [concatsig;svsig];
audiowrite('avqi2.wav',avqiwav2,fs)

if isunix == 1; [s w] = system('Praat.exe --run praat_avqi2.praat avqi2.wav');
elseif ispc == 1; [s w] = system('Praat.exe --run praat_avqi2.praat avqi2.wav');
end
datatemp = importdata('avqi2.txt');
outavqi2.cpps = datatemp(1);
outavqi2.hnr = datatemp(2);
outavqi2.shim = datatemp(3);
outavqi2.shdb = datatemp(4);
outavqi2.slope = datatemp(5);
outavqi2.tilt = datatemp(6);
outavqi2.avqi = datatemp(7);
delete('avqi2.txt')
delete('avqi2.wav')
kjk
function [blah fc fm] = bandlevels(dBspec,f,bands)

% This function calculates the octave-band or third-octave band levels of
% an fft spectrum by adding up frequency bins of certain bands.
% 
% [blah fc fm] = bandlevels(dBspec,f,bands)
%    Inputs: 'dBspec' is the fft dB spectrum of the signal, 'f' is the
%    frequency vector, and 'bands' is the group of bandwidths desired in the output:
%       'iecoct' for all separate IEC octaves (standard octave band analysis)
%       'iec3oct' for all IEC third octaves (standard third-octave band analysis)
%       'iecoctHF' for HFE octaves and summed lower octaves
%       'iec3octHF' for only HFE third octaves
%       'erb' for all separate erbs 1:42
%       'erbg' for erb groupings correlated with octave groupings
%       'oct' for all separate octaves 1:9 (Ternstrom)
%       'octg' for grouped octaves 2:5,6:7,8, and 9
% 
%   Outputs: blah is a 1xN vector containing band values, fc is a 1x(N+1)
%   vector containing the associated cut-on and cut-off frequencies of
%   bands used for the analysis, fm is a 1xN vector containing the
%   associated center (or midband) frequencies for each band if based on
%   the IEC standard. ALL FREQUENCY VALUES ARE GIVEN IN KHZ.
% 
% Example:
%     levels = bandlevels(sigspectrum,sigf,'iecoct');
%     [levels fc fm] = bandlevels(sigspectrum,sigf,'iecoct');

% close all

% create vector with cutoff frequencies for bands of interest

switch bands
    case 'iecoct'
        [junk fm] = octaves(26:34,1,2);
        fc = junk(:,3)*1000; fc(end+1) = 22000;
        fc = fc';
    case 'iecoctHFE'
        [junk fm] = octaves([26 33],1,2);
        fc = junk(:,3)*1000; fc(end+1) = 22000;
        fc = fc';
    case 'iec3oct'
        [junk fm] = octaves(17:43,3,2);
        fc = junk(:,3)*1000; fc(end+1) = 22000;
        fc = fc';
    case 'iec3octHFE'
        [junk fm] = octaves(38:43,3,2);
        fc = junk(:,3)*1000; fc(end+1) = 22000;
        fc = fc';
    case 'erb'
        [junk fc] = erb(1:42);
    case 'erbg'
        [junk fc] = erb([2 17 29 35 42]);
    case 'oct'
        fc = [80 160 320 640 1250 2500 5000 10000 20000];
    case 'octg'
        fc = [80 1250 5000 10000 20000];
end

fci = ceil(fc/f(2));

linspec = 10.^(dBspec/10); % convert to intensity

% add up intensities in separate frequency bands
for ind = 1:length(fci)-1
    if fci(ind) == fci(ind+1)
        fci(ind) = fci(ind) - 1;
    end
    try
        blah(ind) = sum(linspec(fci(ind)+1:fci(ind+1),1:end),1);
    catch
        blah(ind) = sum(linspec(fci(ind)+1:end,1:end),1);
    end
end

% % add up intensities to get total intensity of spectrum
% blah.tot = sum(blah.lin(fci(1)+1:fci(end),1:end),1);

% convert INTENSITIES back to dB (SPL)
blah = 10*log10(blah);
fc = fc/1000; % convert to kHz

% % -------------plotting------------------------
% labels = txt(11,2:end);
% for i=1:length(fci)-1
%     levels(i,:) = [80 75 70 65];
% end
% % figure; plot(blah.f/1000,blah.dB); legend(labels)
% % figure; semilogx(blah.f,blah.dB); legend(labels)
% 
% % calculate dB levels re: to total intensity level
% for ind = 1:length(fci)-1
%     blah.rel(ind,:) = blah.bandsdB(ind,:) - blah.totdB;
% end
% if lev == 1
%     blah.abs = blah.rel + levels;
%     blah.totdBabs = blah.totdB + [80 75 70 65];
% else
%     blah.abs = blah.bandsdB;
% end
% 
% % blah.rel(tot+1,:)=num(10,2:end)-max(num(10,2:end));
% % blah.rel = blah.octdB - overall;

function totalLoudness = specificLoudness2TotalLoudness(specificLoudness,fcAll)
% This function calculates total loudness by integrating specific loudness over ERBrate frequency.
% Input:
%    specificLoudness -- Mx1 real vector, sone/ERBrate, specific loudness function, M auditory
%                        filters
%    fcAll -- Mx1 real vector, Hz, auditory filter center frequencies
% Output:
%    totalLoudness -- real scalar, sone, total loudness of input

% Mark Skowronski, February 12, 2013

% Convert fcAll to ERBrate, calculate delta frequency:
df = diff(f2ERBrate(fcAll)); % ERBrate
df(end+1) = df(end); % repeat last value to match vector lengths

% Integrate:
totalLoudness = sum(specificLoudness.*df);

return;

% Bye!
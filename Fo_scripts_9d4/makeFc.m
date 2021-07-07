function fcAll = makeFc(fcRange,ERBspacing)
% This function creates a vector of center frequencies within fcRange equally spaced in ERB rate
% according to ERBspacing.
% Input:
%    fcRange -- 1x2 real vector, Hz, [min max] frequencies of center frequency range
%    ERBspacing -- real scalar, ERB/ERB, spacing between adjacent center frequencies in ERB rate
% Output:
%    fcAll -- Mx1 real vector, Hz, M center frequencies
% Note: fcAll(end) may undershoot range max because range max in ERB is not integer multiple of
% ERBspacing above range min.

% Convert fcRange to ERBrate:
fcRangeERBrate = f2ERBrate(fcRange); % [min,max] in ERBrate

% Construct fcAll in ERBrate:
fcAllERBrate = [fcRangeERBrate(1):ERBspacing:fcRangeERBrate(2)]'; % COLUMN vector

% Convert from ERBrate to linear frequency:
fcAll = ERBrate2f(fcAllERBrate);

return;

% Bye!
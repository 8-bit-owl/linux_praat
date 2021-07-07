function S = calculateSharpness(specificLoudness,fcAll,gHeight,gBreakpoint)
% This function calculates sharpness according to Zwicker 1999.
% Input [defaults]:
%    specificLoudness -- Mx1 real vector, sones/ERB, specific loudness function
%    fcAll -- Mx1 real vector, Hz, center frequencies of auditory filters
%    gHeight -- real scalar, unitless, height above unity of g weighting function [3]
%    gBreakpoint -- real scalar, barks, breakpoint of g weighting function [16]
% Output:
%    S -- real scalar, aspers, sharpness of specific loudness function

% Mark Skowronski, February 8, 2013

% Check inputs:
if nargin<4
   gBreakpoint = 16; % bark
end;
if nargin<3
   gHeight = 3; 
end;

% Convert fcAll to bark:
z = 13*atan(7.6e-4*fcAll)+3.5*atan((fcAll/7500).^2); % Hz to Bark

% Sharpness weighting function:
g = ones(size(z));
g(z>gBreakpoint) = gHeight/((24-gBreakpoint)^2)*(z(z>gBreakpoint)-gBreakpoint).^2+1;

% Integrate numerator and denominator of sharpness:
dz = diff(z);
dz(end+1) = dz(end); % copy last delta z
num = sum(specificLoudness.*dz.*g.*z);
den = sum(specificLoudness.*dz);

% Calculate sharpness:
S = 0.11*num/den;

return;

% Bye!
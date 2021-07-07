function z = f2bark(f)
% This function converts frequencies f in Hz to critical band rate in Bark.
% Input:
%    f -- Mx1 real vector, Hz, frequencies
% Output:
%    z -- Mx1 real vector, Bark, converted frequencies
%
% Reference: Zwicker and Terhardt, JASA 68(5), pp. 1523-1525, Nov. 1980

% Mark Skowronski, May 20, 2013

%z = 13*atan(0.76*f/1000)+3.5*(atan(f/7500)).^2; % incorrect interpretation of Eq. 1
z = 13*atan(0.76*f/1000)+3.5*atan((f/7500).^2);

return;

% Bye!
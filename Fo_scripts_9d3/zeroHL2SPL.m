function [dBSPL,f] = zeroHL2SPL(f,transducer)
% This function returns the dB SPL level for 0 dB HL at frequencies f according to ANSI S3.6-1996
% for a given transducer.  Cubic spline interpolation is used for frequencies between reference
% frequencies, and cubic spline extrapolation is used for frequencies outside the reference
% frequency range.
% Input: [defaults]
%    f -- Nx1 real vector, Hz, frequencies to return dB SPL level,
%         [125,250,500,750,1000,2000,3000,4000,6000,8000]
%    transducer -- string, 'TDH39', 'TDH49/50', ['ER3A'], 'speaker0bin', 'speaker0mon',
%                  'speaker45mon', 'speaker90mon'
% Output:
%    dBSPL -- Nx1 real vector, dB SPL, value of 0 dB HL at f
%    f -- Nx1 real vector, same as input, or defaults

% Set ANSI frequencies and dB SPL levels:
fRef = [125,250,500,750,1000,1500,2000,3000,4000,6000,8000]'; % Hz, COLUMN vector
dBSPLRef(1).transducer = 'TDH39';
dBSPLRef(1).dBSPL = [45.0,25.5,11.5,8.0,7.0,6.5,9.0,10.0,9.5,15.5,13.0]'; % dB SPL, COLUMN
dBSPLRef(2).transducer = 'TDH49/50';
dBSPLRef(2).dBSPL = [47.5,26.5,13.5,8.5,7.5,7.5,11.0,9.5,10.5,13.5,13.0]'; % dB SPL, COLUMN
dBSPLRef(3).transducer = 'ER3A';
dBSPLRef(3).dBSPL = [26.0,14.0,5.5,2.0,0.0,2.0,3.0,3.5,5.5,2.0,0.0]'; % dB SPL, COLUMN
dBSPLRef(4).transducer = 'speaker0bin';
dBSPLRef(4).dBSPL = [22.0,11.0,4.0,2.0,2.0,0.5,-1.5,-6.0,-6.5,2.5,11.5]'; % dB SPL, COLUMN
dBSPLRef(5).transducer = 'speaker0mon';
dBSPLRef(5).dBSPL = [24.0,13.0,6.0,4.0,4.0,2.5,0.5,-4.0,-4.5,4.5,13.5]'; % dB SPL, COLUMN
dBSPLRef(6).transducer = 'speaker45mon';
dBSPLRef(6).dBSPL = [23.5,12.0,3.0,0.5,0.0,-1.0,-2.5,-9.0,-8.5,-3.0,8.0]'; % dB SPL, COLUMN
dBSPLRef(7).transducer = 'speaker90mon';
dBSPLRef(7).dBSPL = [23.0,11.0,1.5,-1.0,-1.5,-2.5,-1.5,-6.5,-4.0,-5.0,5.5]'; % dB SPL, COLUMN

% Check inputs:
if nargin<2,
   transducer = 'ER3A';
end;
if nargin<1, % return reference frequencies and dBSPL for ER3A
   f = fRef; % Hz
   dBSPL = dBSPLRef(3).dBSPL; % dB SPL, ER3A default
   return;
end;

% Get dB reference for transducer:
t = strmatch(transducer,{dBSPLRef.transducer},'exact');
if isempty(t), % no match, so use default
   dBRef = dBSPLRef(3).dBSPL; % dB SPL, ER3A
else
   dBRef = dBSPLRef(t).dBSPL; % dB SPL
end;

% Use cubic interpolation/extrapolation to determine dBSPL at f:
dBSPL = interp1(fRef,dBRef,f,'spline','extrap');

end


function x = synthCalls01(F0,F1,L,gk,tk,fs),
% function x = synthCalls01(F0,F1,L,gk,tk,fs),
% This function creates piecewise exponential synthetic bat calls
% Input:
%       F0 -- start frequency of call, Hz (default: 30e3)
%       F1 -- end frequency of call, Hz (15e3)
%       L -- call duration, sec (10e-3)
%       gk -- normalized knee frequency, 0<gk<1 (0.5)
%       tk -- normalized time of knee, 0<tk<1 (0.5)
%       fs -- sampling rate, Hz (250e3)
% Output:
%       x -- 1xround(L*fs) vector, time domain of synthetic call
%
% Note: Fknee = Fmin + gk*(Fmax-Fmin) and Tknee = tk*L.  The knee frequency
%       is the frequency of the point on the FM curve furthest from the line
%       between (0,F0) and (L,F1) representing linear frequency modulation.  The
%       knee is also the point on the FM curve with slope equal to the slope of
%       the line representing linear frequency modulation.

% Mark Skowronski, May 17, 2007
% Version 1: Based on synthExp01.m, this file transform the code into function form.
% synthExp01.m: This file performs the synthetic auto-detection experiment.
% Mark Skowronski, May 16, 2007

% Init output:
x = [];

% Check inputs:
if nargin<6,
   fs=250e3; % Hz
end;
if nargin<5,
   tk = .5;
end;
if nargin<4,
   gk = .5;
end;
if nargin<3,
   L = 10e-3; % sec
end;
if nargin<2,
   F1 = 15e3; % Hz
end;
if nargin<1,
   F0 = 30e3; % Hz
end;
if tk<=0 || tk>=1,
   error('ERROR: Must have 0<tk<1.');
   return;
end;
if gk<=0 || gk>=1,
   error('ERROR: Must have 0<gk<1.');
   return;
end;
if fs<=0,
   error('ERROR: Must have fs>0.');
   return;
end;
if L<=0,
   error('ERROR: Must have L>0.');
   return;
end;
if F1<0 || F1>fs/2,
   error('ERROR: Must have 0<=F1<=fs/2.');
   return;
end;
if F0<0 || F0>fs/2,
   error('ERROR: Must have 0<=F0<=fs/2.');
end;

% Determine call type:
%   Type I:   F0<F1 and tk>gk  -- positive-slope concave-up call
%   Type II:  F0>F1 and tk+gk<1  -- negative-slope concave-up call
%   Type III: F0>F1 and tk+gk>1  -- negative-slope concave-down call
%   Type IV:  F0<F1 and tk<gk  -- positive-slope concave-down call
if F0==F1,
   callType = 0; % pure tone (no frequency modulation)
elseif (F0<F1 && gk==tk)|| (F0>F1 && (gk+tk)==1),
   callType = -1; % linear FM
elseif F0<F1 && tk>gk,
   callType = 1; % Type I piecewise exponential
elseif F0>F1 && (tk+gk)<1,
   callType = 2; % Type II piecewise exponential
elseif F0>F1 && (tk+gk)>1,
   callType = 3; % Type III piecewise exponential
else
   callType = 4; % Type IV piecewise exponential
end;

% Create call:
switch callType
   case 0,
      % Pure tone:
      t = [0:round(L*fs)]/fs; % sec
      x = sin(2*pi*F0*t);
   case -1,
      % Linear chirp:
      % f(t) = F0 + (F1-F0)*t/L
      t = [0:round(L*fs)]/fs; % sec
      phi = F0*t+(F1-F0)*t.^2/(2*L);
      x = sin(2*pi*phi);
   case 1,
      % Type I: positive-slope concave-up call
      % Convert gk, tk to universal-type variables:
      gk1 = gk;
      tk1 = tk;
      
      % Solve for piecewise-exponential parameters:
      a1 = fzero(@(a1) exp(a1)*(tk1-gk1*a1)-tk1,[1e-10 400]); % a1>0, always
      a2 = fzero(@(a2) (1-tk1)*(exp(a2)-1)-a2*(1-gk1),[1e-10 400]); % a2>0, always
      
      % f(t) = F0 + (F1-F0)*g(t)
      % g(t) = g1(t), t=[0,tk1*L]
      %      = g2(t), t=[tk1*L,L]
      % g1(t) = gk1*(exp(a1*(t/L)/tk1)-1)/(exp(a1)-1)
      % g2(t) = gk1 + (1-gk1)*(exp(a2*((t/L)-tk1)/(1-tk1))-1)/(exp(a2)-1)
      t1 = [0:round(tk1*L*fs)]/fs; % sec
      phi1 = F0*t1+(F1-F0)*gk1/(exp(a1)-1)*(L*tk1/a1*(exp(a1*t1/(L*tk1))-1)-t1);
      x1 = sin(2*pi*phi1);
      t2 = [round(tk1*L*fs)+1:round(L*fs)]/fs; % sec
      phi2 = (F0+(F1-F0)*(gk1-(1-gk1)/(exp(a2)-1)))*(t2-tk1*L)+(F1-F0)*(1-gk1)/(exp(a2)-1)*L*(1-tk1)/a2*(exp(a2*(t2/L-tk1)/(1-tk1))-1);
      phi2 = phi2+F0*tk1*L+(F1-F0)*gk1/(exp(a1)-1)*(L*tk1/a1*(exp(a1)-1)-tk1*L); % add phi1(tk1*L) for phase continuity
      x2 = sin(2*pi*phi2);
      t = [t1,t2];
      x = [x1,x2];
   case 2,
      % Type II: negative-slope concave-up call
      % Convert gk, tk to universal-type variables:
      gk1 = gk;
      tk1 = 1-tk;
      
      % Solve for piecewise-exponential parameters:
      a1 = fzero(@(a1) exp(a1)*(tk1-gk1*a1)-tk1,[1e-10 400]); % a1>0, always
      a2 = fzero(@(a2) (1-tk1)*(exp(a2)-1)-a2*(1-gk1),[1e-10 400]); % a2>0, always
      
      % f(t) = F0 + (F1-F0)*(1-g(L-t))
      % g(t) = g1(t), t=[0,tk1*L]
      %      = g2(t), t=[tk1*L,L]
      % g1(t) = gk1*(exp(a1*(t/L)/tk1)-1)/(exp(a1)-1)
      % g2(t) = gk1 + (1-gk1)*(exp(a2*((t/L)-tk1)/(1-tk1))-1)/(exp(a2)-1)
      t1 = [0:round(tk*L*fs)]/fs; % sec
      phi1 = F1*t1+(F0-F1)*(gk1-(1-gk1)/(exp(a2)-1))*t1+(F0-F1)*(1-gk1)/(exp(a2)-1)*L*(tk1-1)/a2*(exp(a2*((L-t1)/L-tk1)/(1-tk1))-exp(a2));
      x1 = sin(2*pi*phi1);
      t2 = [round(tk*L*fs)+1:round(L*fs)]/fs; % sec
      phi2 = (F0+(F1-F0)*(1+gk1/(exp(a1)-1)))*(t2-tk*L)+(F1-F0)*(-gk1)/(exp(a1)-1)*(-L*tk1/a1)*(exp(a1*(L-t2)/(tk1*L))-exp(a1*(L-tk*L)/(tk1*L)));
      phi2 = phi2+F1*tk*L+(F0-F1)*(gk1-(1-gk1)/(exp(a2)-1))*tk*L+(F0-F1)*(1-gk1)/(exp(a2)-1)*L*(tk1-1)/a2*(exp(a2*((L-tk*L)/L-tk1)/(1-tk1))-exp(a2));
      x2 = sin(2*pi*phi2);
      t = [t1,t2];
      x = [x1,x2];
   case 3,
      % Type III: negative-slope concave-down call
      % Convert gk, tk to universal-type variables:
      gk1 = 1-gk;
      tk1 = tk;
      
      % Solve for piecewise-exponential parameters:
      a1 = fzero(@(a1) exp(a1)*(tk1-gk1*a1)-tk1,[1e-10 400]); % a1>0, always
      a2 = fzero(@(a2) (1-tk1)*(exp(a2)-1)-a2*(1-gk1),[1e-10 400]); % a2>0, always
      
      % f(t) = F0 + (F1-F0)*g(t)
      % g(t) = g1(t), t=[0,tk1*L]
      %      = g2(t), t=[tk1*L,L]
      % g1(t) = gk1*(exp(a1*(t/L)/tk1)-1)/(exp(a1)-1)
      % g2(t) = gk1 + (1-gk1)*(exp(a2*((t/L)-tk1)/(1-tk1))-1)/(exp(a2)-1)
      t1 = [0:round(tk1*L*fs)]/fs; % sec
      phi1 = F0*t1+(F1-F0)*gk1/(exp(a1)-1)*(L*tk1/a1*(exp(a1*t1/(L*tk1))-1)-t1);
      x1 = sin(2*pi*phi1);
      t2 = [round(tk1*L*fs)+1:round(L*fs)]/fs; % sec
      phi2 = (F0+(F1-F0)*(gk1-(1-gk1)/(exp(a2)-1)))*(t2-tk1*L)+(F1-F0)*(1-gk1)/(exp(a2)-1)*L*(1-tk1)/a2*(exp(a2*(t2/L-tk1)/(1-tk1))-1);
      phi2 = phi2+F0*tk1*L+(F1-F0)*gk1/(exp(a1)-1)*(L*tk1/a1*(exp(a1)-1)-tk1*L); % add phi1(tk1*L) for phase continuity
      x2 = sin(2*pi*phi2);
      t = [t1,t2];
      x = [x1,x2];
   case 4,
      % Type IV: positive-slope concave-down call
      % Convert gk, tk to universal-type variables:
      gk1 = 1-gk;
      tk1 = 1-tk;
      
      % Solve for piecewise-exponential parameters:
      a1 = fzero(@(a1) exp(a1)*(tk1-gk1*a1)-tk1,[1e-10 400]); % a1>0, always
      a2 = fzero(@(a2) (1-tk1)*(exp(a2)-1)-a2*(1-gk1),[1e-10 400]); % a2>0, always
      
      % f(t) = F0 + (F1-F0)*(1-g(L-t))
      % g(t) = g1(t), t=[0,tk1*L]
      %      = g2(t), t=[tk1*L,L]
      % g1(t) = gk1*(exp(a1*(t/L)/tk1)-1)/(exp(a1)-1)
      % g2(t) = gk1 + (1-gk1)*(exp(a2*((t/L)-tk1)/(1-tk1))-1)/(exp(a2)-1)
      t1 = [0:round(tk*L*fs)]/fs; % sec
      phi1 = F1*t1+(F0-F1)*(gk1-(1-gk1)/(exp(a2)-1))*t1+(F0-F1)*(1-gk1)/(exp(a2)-1)*L*(tk1-1)/a2*(exp(a2*((L-t1)/L-tk1)/(1-tk1))-exp(a2));
      x1 = sin(2*pi*phi1);
      t2 = [round(tk*L*fs)+1:round(L*fs)]/fs; % sec
      phi2 = (F0+(F1-F0)*(1+gk1/(exp(a1)-1)))*(t2-tk*L)+(F1-F0)*(-gk1)/(exp(a1)-1)*(-L*tk1/a1)*(exp(a1*(L-t2)/(tk1*L))-exp(a1*(L-tk*L)/(tk1*L)));
      phi2 = phi2+F1*tk*L+(F0-F1)*(gk1-(1-gk1)/(exp(a2)-1))*tk*L+(F0-F1)*(1-gk1)/(exp(a2)-1)*L*(tk1-1)/a2*(exp(a2*((L-tk*L)/L-tk1)/(1-tk1))-exp(a2));
      x2 = sin(2*pi*phi2);
      t = [t1,t2];
      x = [x1,x2];
end;
      


% Bye!
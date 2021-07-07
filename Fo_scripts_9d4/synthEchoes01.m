function yCall = synthEchoes01(xCall,micLoc,sourceLoc,parametersEcho,hAAFilters)
% function yCall = synthEchoes01(xCall,micLoc,sourceLoc,parametersEcho,hAAFilters)
% Based on echoFilter01.m, this function adds echoes to a synthetic call.
% Mark Skowronski, October 9, 2007
% echoFilter01.m: This file creates diffuse echoes scattered off the ground from
% a synthetic echolocation call.
% Mark Skowronski, October 6, 2007

if nargin<5,
   % Design hAA filters:
   pathLength = [.1:.1:40]; % m
   hAAFilters = cell(1,length(pathLength));
   filterWindow = hamming(parametersEcho.NaaFilter+1);
   f = [0 [10e3:10e3:120e3]/125e3 1];
   for p=1:length(hAAFilters),
      fGain = -pathLength(p)*[0  0.1594 0.5258 0.9387 1.3204 1.6631 1.9817 2.2923 2.6073 2.9351 3.2813 3.6500 4.0435 4.5]; % dB
      hAAFilters{p} = fir2(parametersEcho.NaaFilter,f,10.^(fGain/20),filterWindow); % ROW vector
   end;
end;
      
% Generate random points on lower-half sphere:
z=-rand(parametersEcho.M,1);
t=rand(size(z))*2*pi;
r=sqrt(1-z.^2);
x=r.*sin(t);
y=r.*cos(t);
randRay = [x,y,z]; % unit length
randRaySpec = [x,y,-z]; % reflected vector

% Find angle between random ray specular reflection and diffuse reflection through micLoc:
diffuseLoc = repmat(sourceLoc,size(z,1),1)+randRay*sourceLoc(3)./repmat(-z,1,3);
vecDiffMic = repmat(micLoc,size(z,1),1)-diffuseLoc;
vecDiffMic = vecDiffMic./repmat(sqrt(sum(vecDiffMic.^2,2)),1,3); % norm length
diffuseAngle = acos(sum(randRaySpec.*vecDiffMic,2)); % acos of inner product
pathLength = sqrt(sum((repmat(sourceLoc,size(z,1),1)-diffuseLoc).^2,2))+...
   sqrt(sum((repmat(micLoc,size(z,1),1)-diffuseLoc).^2,2)); % m, COLUMN vector

% For atmospheric absorption, assume 20 C, 50% RH, 1 atm:
% dB/m: 0.1594 0.5258 0.9387 1.3204 1.6631 1.9817 2.2923 2.6073 2.9351 3.2813 3.6500 4.0435
% kHz:  10     20     30     40     50     60     70     80     90     100    110    120
% From atmLoss01.m.

% Create filter for direct path between source and mic:
f = [0 [10e3:10e3:120e3]/125e3 1];
Lsm = norm(micLoc-sourceLoc); % m
fGain = -Lsm*[0  0.1594 0.5258 0.9387 1.3204 1.6631 1.9817 2.2923 2.6073 2.9351 3.2813 3.6500 4.0435 4.5]; % dB
hAA = fir2(parametersEcho.NaaFilter,f,10.^(fGain/20)); % ROW vector
delayCont = Lsm/parametersEcho.c*parametersEcho.fs; % samples, continuous value
delayInt0 = round(delayCont); % samples, integer value delay
%delayFrac0 = delayCont - delayInt0; % fraction of sample, in range [-.5,+.5], delayInt+delayFrac=delayCont
%sinArgRange = [-60:60]; % samples
%hFD = sin(pi*(sinArgRange-delayFrac0))./(pi*(sinArgRange-delayFrac0)); % ROW vector
%hC = conv(hAA,hFD); % ROW vector
hC = hAA; % ROW vector, ignore fractional delay, for program speed
hS0 = hC*parametersEcho.L0/Lsm;

% For each diffuse reflection close to mic, design filter:
delayCont = pathLength/parametersEcho.c*parametersEcho.fs; % samples, continuous value
delayInt = round(delayCont); % samples, integer value delay
%delayFrac = delayCont - delayInt; % fraction of sample, in range [-.5,+.5], delayInt+delayFrac=delayCont
pIndex = find(diffuseAngle<parametersEcho.sigDif);
%sinArg = pi*(sinArgRange(ones(length(pIndex),1),:)-delayFrac(pIndex,ones(1,length(sinArgRange))));
%hFD = sin(sinArg)./sinArg; % ROW vectors
hCScaleFactor = exp(-diffuseAngle(pIndex).^2/(.5*parametersEcho.sigDif^2))*...
   parametersEcho.L0./pathLength(pIndex)*10^(-parametersEcho.groundAbs/20); % off-angle + spherical spread + ground loss
hS = cell(1,length(pIndex)); % scaled combination of AA and FD filters
for p=1:length(pIndex),
   % Get hAA filter:
   hAA = hAAFilters{ceil(pathLength(pIndex(p))*10)};

   % Combine filters, scale:
%   hC = conv(hAA,hFD(p,:)); % ROW vector
   hC = hAA; % ROW vector, ignore fractional delay, for program speed
   hS{p} = hC*hCScaleFactor(p);
end;

% Combine hS0 with hS using proper delays:
%hTotal = zeros(1,length(hS0)+max(delayInt(diffuseAngle<parametersEcho.sigDif)));
hTotal = zeros(1,length(hS0)+max(delayInt(diffuseAngle<parametersEcho.sigDif))-delayInt0); % ignore leading zeros
hTotalRange = [1:length(hS0)];
%hTotal(hTotalRange+delayInt0) = hTotal(hTotalRange+delayInt0)+hS0;
hTotal(hTotalRange+0) = hTotal(hTotalRange+0)+hS0;
for p=1:length(pIndex),
%   hTotal(hTotalRange+delayInt(pIndex(p))) = hTotal(hTotalRange+delayInt(pIndex(p)))+hS{p};
   hTotal(hTotalRange+delayInt(pIndex(p))-delayInt0) = hTotal(hTotalRange+delayInt(pIndex(p))-delayInt0)+hS{p};
end;

% Filter:
yCall = conv(hTotal,xCall);
yCall = [zeros(1,delayInt0),yCall]; % append leading zeros (propagation delay from source to mic)

% Bye!
function x = xFilter(x,z,p,k)
% This function inputs poles/zeros/gain term for a filter, constructs 1st- and 2nd-order systems,
% and filters x using filterSys.

% Convert poles/zeros/gain term into S1 and S2:
[S1,S2] = getS1S2(z,p);

% Filter:
x = filterSys(S1,S2,x);

% Scale by k:
x = x*k;

return;



function [S1,S2] = getS1S2(z,p)
% This function converts the poles/zeros to S1/S2 matrices.

% Find real/complex poles/zeros (may be empty):
[pReal,pC] = getRealCC(p);
[zReal,zC] = getRealCC(z);

% Combine complex poles/zeros with their CCs (may be empty):
b012 = [ones(length(zC),1),-2*real(zC),real(zC).^2+imag(zC).^2];
a12 = [-2*real(pC),real(pC).^2+imag(pC).^2];

% Combine real poles/zeros to form 2nd-order terms:
if mod(length(z),2)==1, % odd length: 1 S1 filter, rest in S2
   % Construct S1, make index vector:
   S1 = [1,-zReal(1),-pReal(1)];
   tb = [2:2:length(zReal)];
   ta = [2:2:length(pReal)];
else % even length: no S1, all in S2
   S1 = [];
   tb = [1:2:length(zReal)];
   ta = [1:2:length(pReal)];
end;
   
% Construct S2 from b012, a12, and remaining real poles/zeros:
b012 = [b012;[ones(length(tb),1),-(zReal(tb)+zReal(tb+1)),zReal(tb).*zReal(tb+1)]];
a12 = [a12;[-(pReal(ta)+pReal(ta+1)),pReal(ta).*pReal(ta+1)]];
S2 = [b012,a12];

return;



function [xReal,xC] = getRealCC(x)

% Get real terms:
xReal = x(imag(x)==0); % may be empty

% Get all complex terms with positive imaginary part (assume CC exists):
xC = x(imag(x)>0); % may be empty

return;

function y = filterSys(S1,S2,x)
% This function filters x with the first- and second-order filters in S1, S2, respectively.
% Inputs:
%  S1 -- N1x3 real matrix, filter coefficients for all N1 first-order filters: [b0 b1 a1] each row
%  S2 -- N2x5 real matrix, filter coefficients for all N2 second-order filters: [b0 b1 b2 a1 a2]
%  x -- T-length vector, real values
% Output:
%  y -- T-length vector, real values
%
% First-order filters:
%  y(n) = b0*x(n) + b1*x(n-1) - a1*y(n-1)
%       = S1(i,1)*x(n) + S1(i,2)*x(n-1) - S1(i,3)*y(n-1), i^th filter
% Second-order filters:
%  y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2)
%       = S2(i,1)*x(n) + S2(i,2)*x(n-1) + S2(i,3)*x(n-2) - S2(i,4)*y(n-1) - S2(i,5)*y(n-2)

% Mark Skowronski, August 10, 2009

% Init output:
y = x;

% Check inputs:
if nargin<3,
   warning('WARNING: requires S1, S2, and x input.');
   return;
end;
if isempty(x), % nothing to filter, so return empty output
   return;
end;
if isempty(S1) && isempty(S2), % nothing to filter with, so return input
   return;
end;
if ~isempty(S2),
   if size(S2,2)~=5,
      warning('WARNING: S2 must have 5 columns if not empty.');
      return;
   end;
end;
if ~isempty(S1),
   if size(S1,2)~=3,
      warning('WARNING: S1 must have 3 columns if not empty.');
      return;
   end;
end;

% Filter with each S1 filter:
for pN1 = 1:size(S1,1),
   b = S1(pN1,1:2);
   a = [1,S1(pN1,3)]; % include a0
   y = filter(b,a,y); % replace with filter1() when porting
end;

% Filter with each S2 filter:
for pN2 = 1:size(S2,1),
   b = S2(pN2,1:3);
   a = [1,S2(pN2,4:5)]; % include a0
   y = filter(b,a,y); % replace with filter1() when porting
end;

return;




 % Bye!
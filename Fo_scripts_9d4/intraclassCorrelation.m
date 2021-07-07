function [r,rApprox] = intraclassCorrelation(x)
% This function calculates the intraclass correlation of x.
% Input:
%    x -- NxJ real matrix, N data points, J measures (ratings) per data point
% Output:
%    r -- intraclass correlation of x
%    rApprox -- large-J approximation of r
% Note: For ratings, the J measures may be 1 measure from J raters (inter-rater) or J measures from
% 1 rater (intra-rater).

% Mark Skowronski, March 7, 2013

% Get size of x:
[N,J] = size(x);

% Calculate mean of all x values:
xBar = mean(x(:));

% Remove mean from x:
x0 = x-xBar;

% Calculate unbiased variance of x:
s2 = var(x0(:));

% Calculate intraclass correlation:
r = 0; % init
for n=1:N,
   for j=1:J-1,
      for i=j+1:J,
         r = r + x0(n,i)*x0(n,j);
      end;
   end;
end;

% Normalize r:
df = N*J*(J-1)/2-1;
r = r/(df*s2);

% Alternative r value (large J approximation):
xn = mean(x,2); % mean along J measures
rApprox = J/(J-1)*sum((xn-xBar).^2)/(N*s2)-1/(J-1);

return;

% Bye!
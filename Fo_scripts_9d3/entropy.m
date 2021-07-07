function H = entropy(f)
% Calculate entropy of a discrete distribution
% Usage: H = entropy(f)
%  f - input distribution as a vector
%  H - entropy
N = length(f);
H = 0;
for j = 1:N
   H = H - f(j) * logz(f(j));
end

function y = logz(x)
if (x > 0)
   y = log(x);
else
   y = 0;
end

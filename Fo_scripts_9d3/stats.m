function s = stats(x,trim)
% Statistical measures of a vector. The input vector should not contain any
% zeros unless they are truly part of the data. This will never be the case
% for a frequency vector. Should add error checking at some point.
% x = round(rand(1,100)*100);
% trim = 5;

% The trim is a percentage of elements removed from each end of the sorted
% array x. If trim = 10% and length(x) = 100, then 10 elements are removed
% from each end corresponding to the 10 lowest/highest.
if exist('trim')
    if trim < 0 trim = abs(trim); end
    if trim > 1 trim = trim/100; end
    t = round(trim*length(x));
end

% Pre-allocate the stats structure
s = struct('n',[],'min',[],'max',[],'mean',[],'var',[],'std',[],...
    'median',[],'Q1',[],'Q3',[],'IQR',[],'mode',[],'skew',[],...
    'bowley_skew',[],'trim',[],'trim_n',[],'trim_min',[],'trim_max',[],...
    'trim_mean',[],'trim_var',[],'trim_std',[]);

% fill the stats structure
s.n = length(x);
s.min = min(x);
s.max = max(x);
s.mean = mean(x);
s.var = var(x);
s.std = std(x);
s.median = median(x);
s.Q1 = median(x(find(x <= s.median)));
s.Q3 = median(x(find(x >= s.median)));
s.IQR = s.Q3 - s.Q1;
s.mode = mode(x);
s.skew = (s.n/((s.n-1)*(s.n-2)*(s.std^3)))*sum((x-s.mean).^3);
if s.IQR ~= 0
    s.bowley_skew = (s.Q1-2*s.median+s.Q3)/s.IQR;
else
    s.bowley_skew = NaN;
end

% trimmed stats - useful for data with extreme outliers
if (exist('t') && t ~= 0)
    xs = sort(x);
    xt = xs(t+1:end-t);
    s.trim = trim;
    s.trim_n = length(xt);
    s.trim_min = min(xt);
    s.trim_max = max(xt);
    s.trim_mean = median(xt);
    s.trim_var = var(xt);
    s.trim_std = std(xt);
end


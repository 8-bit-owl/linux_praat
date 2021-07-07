function s=basicstats(x)
% basicstats
% this function returns the basic stats of an array in a structure called
% 'stats'.  some stuff from Albert's script and some from
%     D.C. Hanselman, University of Maine, Orono, ME 04469
%     MasteringMatlab@yahoo.com
%
% S.n = number of observations (excluding NaNs)
% S.nan = number of NaNs
% S.min = minimum
% S.q1 = 25th percentile
% S.median = median (50th percentile)
% S.q3 = 75th percentile
% S.max = maximum
% S.mode = mode
% S.mean = mean
% S.std = standard deviation
% S.var = variance
% S.skew = skewness
% S.kurt = kurtosis (= 3 for Gaussian or Normal Data)
% S.p01 = first optional percentile
% S.p02 = second optional percentile, etc.
%
%from http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
% "Skewness is a measure of symmetry, or more precisely, the lack of 
% symmetry. A distribution, or data set, is symmetric if it looks the same
% to the left and right of the center point. ... The skewness for a normal
% distribution is zero, and any symmetric data should have a skewness 
% near zero. Negative values for the skewness indicate data that are 
% skewed left and positive values for the skewness indicate data that are
% skewed right. By skewed left, we mean that the left tail is long 
% relative to the right tail. 
% 
% "Kurtosis is a measure of whether the data are peaked or flat relative
% to a normal distribution. That is, data sets with high kurtosis tend to
% have a distinct peak near the mean, decline rather rapidly, and have 
% heavy tails. Data sets with low kurtosis tend to have a flat top near the
% mean rather than a sharp peak. A uniform distribution would be the
% extreme case.  ... This definition is used so that the standard normal
% distribution has a kurtosis of zero. In addition, with the second 
% definition positive kurtosis indicates a "peaked" distribution and
% negative kurtosis indicates a "flat" distribution."
%
% eric hunter 20071017

s = struct('n',[],'min',[],'max',[],'mean',[],'var',[],'std',[],...
  'median',[],'Q1',[],'Q3',[],'IQR',[],'mode',[],'skew',[],...
  'bowley_skew',[],'kurt',[]);

if ~isempty(x)

  % fill the stats structure
  s.n =  length(x);
  s.min =  min(x);
  s.max =  max(x);
  s.mean =  sum(x)/s.n;
%   s.mean =  (mean(x));

  %% one way
  xs=x-s.mean;
  xss=xs.*xs;
%   s.var=mean(xss);                                        % var
  s.var=sum(xss)/s.n;                                        % var
  s.std=sqrt(s.var);                                      % std
%   s.skew=mean(xss.*xs)./s.var.^(1.5);                     % skew
%   s.kurt=mean(xss.*xss)./s.var.^2;                        % kurtosis
  s.skew=(sum(xss.*xs)/s.n)./s.var.^(1.5);                % skew
  s.kurt=(sum(xss.*xss)/s.n)./s.var.^2;                   % kurtosis

  
  
%   %% the other way
%   s.var =  (var(x));
%   s.std =  (std(x));
%   s.skew =  ((s.n/((s.n-1)*(s.n-2)*(s.std^3)))*sum((x-s.mean).^3));

  s.median =  (median(x));
  s.Q1 =  (median(x(x <= s.median)));
  s.Q3 =  (median(x(x >= s.median)));
  s.IQR =  (s.Q3 - s.Q1);
  s.mode =  (mode(round(x)));
  
  if s.IQR ~= 0
    s.bowley_skew =  ((s.Q1-2*s.median+s.Q3)/s.IQR);
  else
    s.bowley_skew = NaN;
  end

% xs=sort(x);         % sort data down columns, this puts NaNs last
% 
% tmp=nan+zeros(size(s.mean),class(x));   % template
% s.q1_b=tmp;
% s.median=tmp;
% s.q3_b=tmp;
  
else

  s.n =    0;       %((0));
  s.min =  0;       %(min(0));
  s.max =  0;       %(max(0));
  s.mean = 0;       %(mean(0));
  s.var =  0;       %(var(0));
  s.std =  0;       %(std(0));
  s.median = 0;     % (median(0));
  s.Q1 =   0;
  s.Q3 =   0;
  s.IQR =  0;       %(s.Q3 - s.Q1);
  s.mode = 0;       % (mode(0));
  s.skew = NaN;
  s.bowley_skew = NaN;
  s.kurt=0;

end
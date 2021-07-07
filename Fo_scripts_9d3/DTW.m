function [g,minPath] = DTW(c1,c2,showPlot)
% This function performs dynamic time warping between matrices c1 and c2.  A Type I local constraint
% is used. Euclidean distance is used between feature vectors in c1, c2. See Rabiner and Juang 1993
% for details.
% Input:
%    c1 -- LxN1 real matrix, L features, N1 frames
%    c2 -- LxN2 real matrix, L features, N2 frames
%    showPlot -- 1==show distance matrix w/ minPath; 0==don't (default)
% Output:
%    g -- real scalar, minimum path distance
%    minPath -- Nx2 integer matrix, [i1,i2] pair of frame indices into c1 and c2 (respectively),
%               N = max(N1,N2)

% Mark Skowronski, April 19, 2012

% Check inputs:
if size(c1,1)~=size(c2,1),
   error('ERROR: c1, c2 must have same number of ROWS.'); % Returns w/o continuing
end;
if nargin < 3,
   showPlot = 0; % Default, don't show distance matrix w/ minPath
end;

% Initialize program parameters:
slopeFactor = 2; % Controls the exclusion regions of the global search space.  Larger values exclude less.

% Find Euclidean distance matrix between c1 and c2:
N1 = size(c1,2);
N2 = size(c2,2);
d = zeros(N1,N2); % distance matrix
for p=1:N1, % for each frame of c1
   % Get frame from c1:
   x = c1(:,p);
   
   % Find Euclidean distance between x and all frames in c2:
   d(p,:) = sqrt(sum((c2-repmat(x,1,N2)).^2));
end;
   
% Find start and end columns for each row of search space:
ma = (N1-1)/(N2-1)*slopeFactor; % Slopes
mb = (N1-1)/(N2-1)/slopeFactor;
rowVector = [1:N1];
N2Start = ceil(max((rowVector-1)/ma+1,(rowVector-1-(N1-1))/mb+N2));
N2End = floor(min((rowVector-1)/mb+1,(rowVector-1-(N1-1))/ma+N2));

% Initialize path distance to inf, based on size of d, using local path constraints I:
DI = ones(N1,N2)*inf;
DI(1,1) = d(1,1);

% Construct minimum path matrix DI:
for ir = 2:N1,
   for ic = N2Start(ir):N2End(ir),
      % Use Type I constraint:
      d1 = DI(max(1,ir-1),ic) + d(ir,ic); % at least row 1
      d2 = DI(max(1,ir-1),max(1,ic-1)) + 2*d(ir,ic); % at least row 1, column 1
      d3 = DI(ir,max(1,ic-1)) + d(ir,ic); % at least column 1
      DI(ir,ic) = min([d1,d2,d3]);
   end;
end;

% Save minimum path length:
g = DI(N1,N2);

% Find minimum path through DI:
minPath(1,1:2) = [N1,N2];
for p = 2:N1+N2,
   % Get row/column of previous path element:
   ir = minPath(p-1,1);
   ic = minPath(p-1,2);
   
   % Find out path element that led to ir/ic:
   d1 = DI(max(1,ir-1),ic) + d(ir,ic); % at least row 1
   d2 = DI(max(1,ir-1),max(1,ic-1)) + 2*d(ir,ic); % at least row 1, column 1
   d3 = DI(ir,max(1,ic-1)) + d(ir,ic); % at least column 1
   [junk,a] = min([d1,d2,d3]);
   
   % Update minPath:
   switch a
      case 1
         minPath(p,:) = [max(1,ir-1),ic]; % at least row 1
      case 2
         minPath(p,:) = [max(1,ir-1),max(1,ic-1)]; % at least row 1, column 1
      otherwise
         minPath(p,:) = [ir,max(1,ic-1)]; % at least column 1
   end;
   
   if minPath(p,1)==1 && minPath(p,2)==1,
      break; % Leave for loop
   end;
end;

% Flip, start w/ (1,1) and end w/ (N1,N2):
minPath = flipud(minPath);

% Plot, if desired:
if showPlot,
   figure;
   imagesc(d);
   axis xy;
   hold on;
   plot(minPath(:,2),minPath(:,1),'w-x');
   xlabel('c2');
   ylabel('c1');
end;
   


% Bye!

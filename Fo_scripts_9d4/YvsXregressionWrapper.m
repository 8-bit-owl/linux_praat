function [output,outCell,ANCOVAstats,ANCOVAcoeffs] = YvsXregressionWrapper(x,y,talker,condition,depVariableName)
% This function inputs data and calls YvsXregression() for each talker.
% Input:
%    x -- Nx1 real vector, ind. variable, N samples
%    y -- Nx1 real vector, dep. variable, N samples
%    talker -- Nx1 cell array, talker string
%    condition -- Nx1 cell array, condition string (control, pre, post)
%    depVariableName -- string, name of dependent variable name (used in figure y axis label)
% Output:
%    output -- struct of output info, one element per test talker
%    outCell -- cell array, output info in table form, top row is header, one row per test talker
%    ANCOVAstats -- cell array, ANCOVA results

% Mark Skowronski, March 4, 2014

close all hidden;

% Check inputs:
if nargin<5,
   depVariableName = 'Measure';
end;

% Get training data:
tTrain = strcmp('Control',condition); % logical index
xTrain = x(tTrain);
yTrain = y(tTrain);

% Get non-control talkers:
testTalker = unique(talker(~tTrain));

% Process each non-control talker:
for p=1:length(testTalker),
   % Get pre, post data indices:
   tPre = strcmp('PD Pre-Dx',condition) & strcmp(testTalker{p},talker); % logical index
   tPost = strcmp('PD Post-Dx',condition) & strcmp(testTalker{p},talker); % logical index
   
   % Get pre, post data:
   xPre = x(tPre);
   yPre = y(tPre);
   xPost = x(tPost);
   yPost = y(tPost);
   
   % Regression, store results:
   [output(p),outputCell] = YvsXregression(xTrain,yTrain,xPre,yPre,xPost,yPost,testTalker{p},depVariableName);
   if p==1, % if first talker
      outCell = outputCell; % copy header
   else
      outCell = [outCell;outputCell(2,:)]; % skip header
   end;
end;

% ANCOVA on test population:
[~,ANCOVAstats,ANCOVAcoeffs,statsPop] = aoctool(x,y,condition,.05,'','','','on','separate lines');
figure;
multcompare(statsPop,'estimate','slope');

return;

% Bye!
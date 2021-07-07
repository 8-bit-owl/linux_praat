function [output,outputCell] = YvsXregression(xTrain,yTrain,xPre,yPre,xPost,yPost,talker,depVariableName)
% This function trains a linear regression model then measures the error between model(x) and y for
% two test data conditions, pre and post.
% Input:
%    xTrain,yTrain -- Nx1 real vector, ind. and dep. variables, training data, N samples
%    xPre,yPre -- Mx1 real vector, ind. and dep. variables, pre condition, M samples
%    xPost,yPost -- Lx1 real vector, ind. and dep. variables, post condition, L samples
%    talker -- string, test talker token
%    depVariableName -- string, name of dependent variable measure
% Output:
%    output -- struct of output info
%    outputCell -- cell array of output info, includes header row

% Mark Skowronski, March 4, 2014

% Parameters:
fontName = 'arial';
fontSize = 18;
markerSize = 5;
markerLineWidth = 3;

% Check inputs:
if nargin<8
   depVariableName = 'Measure';
end;
if nargin<7
   talker = 'Talker';
end;

% Get regression of training data:
% Note: b = [a0 a1], stats = [R^2,F,p,error variance].
[bTrain,~,~,~,statsTrain] = regress(yTrain,[ones(length(xTrain),1),xTrain]);

% Get error for each test condition:
errorPre = getError(bTrain,xPre,yPre);
errorPost = getError(bTrain,xPost,yPost);

% Perform t-test for each test condition:
[~,pPre,~,statsPre] = ttest(errorPre);
[~,pPost,~,statsPost] = ttest(errorPost);

% Perform linear regression on error vs. x:
[bPre,~,~,~,statsPre2] = regress(errorPre,[ones(length(xPre),1),xPre]);
[bPost,~,~,~,statsPost2] = regress(errorPost,[ones(length(xPost),1),xPost]);

% Perform t-test comparing each test condition:
[~,pComp,~,statsComp] = ttest2(errorPre,errorPost);

% Perform ANCOVA comparing error for each test condition:
if ~isempty(xPre) && ~isempty(xPost), % some pre and post data
   X = [xPre;xPost];
   Y = [errorPre;errorPost];
   G = [ones(length(xPre),1)*1;ones(length(xPost),1)*2]; % 1=pre, 2=post
   [~,aANC,cANC,statsANC] = aoctool(X,Y,G,.05,'','','','off','separate lines');
else
   aANC = cell(6,6);
   cANC = [];
   statsANC = [];
end;

% Plot results:
figure;
hold on;
plot(xTrain,yTrain,'k.','markersize',10);
xlabel('Age');
ylabel(depVariableName);
title(talker);
plot(xPre,yPre,'go','markersize',markerSize,'linewidth',markerLineWidth);
plot(xPost,yPost,'ro','markersize',markerSize,'linewidth',markerLineWidth);
xLim = [min(xTrain),max(xTrain)];
plot(xLim,bTrain(1)+bTrain(2)*xLim,'k-','linewidth',3);
legend({'Training','Pre','Post'});
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel'),get(gca,'title')],'fontname',fontName,'fontsize',fontSize);

% Store output:
output(1).talker = talker;
output(1).depVariableName = depVariableName;
output(1).regressTrain(1).b = bTrain;
output(1).regressTrain(1).rSquared = statsTrain(1);
output(1).regressTrain(1).F = statsTrain(2);
output(1).regressTrain(1).p = statsTrain(3);
output(1).regressTrain(1).var = statsTrain(4);
output(1).regressTrain(1).N = length(xTrain);
output(1).ttestErrorPre(1).p = pPre;
output(1).ttestErrorPre(1).stats = statsPre;
output(1).ttestErrorPost(1).p = pPost;
output(1).ttestErrorPost(1).stats = statsPost;
output(1).regressErrorPre(1).b = bPre;
output(1).regressErrorPre(1).rSquared = statsPre2(1);
output(1).regressErrorPre(1).F = statsPre2(2);
output(1).regressErrorPre(1).p = statsPre2(3);
output(1).regressErrorPre(1).var = statsPre2(4);
output(1).regressErrorPre(1).N = length(xPre);
output(1).regressErrorPost(1).b = bPost;
output(1).regressErrorPost(1).rSquared = statsPost2(1);
output(1).regressErrorPost(1).F = statsPost2(2);
output(1).regressErrorPost(1).p = statsPost2(3);
output(1).regressErrorPost(1).var = statsPost2(4);
output(1).regressErrorPost(1).N = length(xPost);
output(1).ttestErrorCompare(1).p = pComp;
output(1).ttestErrorCompare(1).stats = statsComp;
output(1).ANCOVA(1).results = aANC;
output(1).ANCOVA(1).coeffs = cANC;
output(1).ANCOVA(1).stats = statsANC;

% Create outputCell from output:
outputCell = makeOutputCell(output);

return;

function errorTerm = getError(bTrain,x,y)
% This function returns the error variable for model(x)-y(x).

% Error = model - y:
errorTerm = (bTrain(1)+bTrain(2)*x) - y;

return;

function outputCell = makeOutputCell(output)
% This function stores output info in cell array form.  Header in top ROW.

outputCell = {...
   'Talker',output.talker;...
   'Dependent variable',output.depVariableName;...
   'Train y-int',output.regressTrain.b(1);...
   'Train slope',output.regressTrain.b(2);...
   'Train r^2',output.regressTrain.rSquared;...
   'Train F',output.regressTrain.F;...
   'Train p',output.regressTrain.p;...
   'Train N',output.regressTrain.N;...
   't-test error Pre d.f.',output.ttestErrorPre.stats.df;...
   't-test error Pre t',output.ttestErrorPre.stats.tstat;...
   't-test error Pre p',output.ttestErrorPre.p;...
   't-test error Post d.f.',output.ttestErrorPost.stats.df;...
   't-test error Post t',output.ttestErrorPost.stats.tstat;...
   't-test error Post p',output.ttestErrorPost.p;...
   't-test error Pre-Post d.f.',output.ttestErrorCompare.stats.df;...
   't-test error Pre-Post t',output.ttestErrorCompare.stats.tstat;...
   't-test error Pre-Post p',output.ttestErrorCompare.p;...
   'Regress error Pre y-int',output.regressErrorPre.b(1);...
   'Regress error Pre slope',output.regressErrorPre.b(2);...
   'Regress error Pre r^2',output.regressErrorPre.rSquared;...
   'Regress error Pre F',output.regressErrorPre.F;...
   'Regress error Pre p',output.regressErrorPre.p;...
   'Regress error Pre N',output.regressErrorPre.N;...
   'Regress error Post y-int',output.regressErrorPost.b(1);...
   'Regress error Post slope',output.regressErrorPost.b(2);...
   'Regress error Post r^2',output.regressErrorPost.rSquared;...
   'Regress error Post F',output.regressErrorPost.F;...
   'Regress error Post p',output.regressErrorPost.p;...
   'Regress error Post N',output.regressErrorPost.N;...
   'ANCOVA group d.f.',output.ANCOVA.results{2,2};...
   'ANCOVA group F',output.ANCOVA.results{2,5};...
   'ANCOVA group p',output.ANCOVA.results{2,6};...
   'ANCOVA ind. var. d.f.',output.ANCOVA.results{3,2};...
   'ANCOVA ind. var. F',output.ANCOVA.results{3,5};...
   'ANCOVA ind. var. p',output.ANCOVA.results{3,6};...
   'ANCOVA interaction d.f.',output.ANCOVA.results{4,2};...
   'ANCOVA interaction F',output.ANCOVA.results{4,5};...
   'ANCOVA interaction p',output.ANCOVA.results{4,6};...
   }'; % header in top ROW

return;

% Bye!
function [output,modelsAll,modelsLOO,cmAll,cmLOO,cmAllSample,cmLOOSample] = crossValidateTestTrain(data,parameters)
% This function performs a cross-validation classification experiment on data.  Each class is
% represented by a Gaussian mixture model.
% Input:
%    data -- Mx1 struct of input data, M data samples
%       .class -- real scalar or string, class identification
%       .talker -- real scalar or string, talker identification
%       .features -- real scalar, ROW vector, or matrix of ROW vectors, measured features for current
%                   data sample
%    parameters -- struct of modeling parameters [defaults]
%       .numMixtures -- integer scalar, number of GMM mixtures for each model [1]
%       .numIterations -- integer scalar, number of GMM training iterations [100]
%       .tol -- real scalar, stopping criterion tolerance [1e-6]
%       .covType -- string, ['full'],'diag' covariance matrix constraint
%       .sharedCov -- logical scalar, true=same covariance matrix for all classes, [false]=unique
%                     covariance matrix for each class
% Output:
%    output -- struct of classification experiment results

% Mark Skowronski, June 10, 2013

% Set defaults:
defaults(1).numMixtures = 1;
defaults(1).numIterations = 100;
defaults(1).tol = 1e-6;
defaults(1).covType = 'full';
defaults(1).sharedCov = false;

% Check inputs:
if nargin<2
   parameters = defaults;
else
   parameters = checkParameters(parameters,defaults);
end;

% Create GMM options, store parameter values:
GMMoptions = statset('gmdistribution'); % default struct
GMMoptions.MaxIter = parameters.numIterations;
GMMoptions.TolFun = parameters.tol;
GMMoptions.Display = 'off';

% Get list of unique classes and talkers:
[classID,talkerID,data] = getDataList(data); % cell arrays of strings

% Train models for each class using ALL data:
for p=1:length(classID), % for each class
   % Get index of training samples:
   k = find(strcmp(classID{p},{data.class})); % index into data
   
   % Pool data:
   x = [data(k(1)).features]; % scalar, ROW vector, or matrix of ROW vectors, grow dynamically
   for p1=2:length(k),
      x = [x;data(k(p1)).features];
   end;
   
   % Train model:
   modelsAll(p).obj = gmdistribution.fit(x,parameters.numMixtures,'covtype',parameters.covType,...
      'sharedcov',parameters.sharedCov,'options',GMMoptions);
   modelsAll(p).class = classID{p};
end;

% Train leave-one-out models that exclude each talker:
for p=1:length(talkerID),
   % Get indices of training samples:
   h1 = strcmp(talkerID{p},{data.talker}); % logical output, case sensitive, indices for current talker
   k1 = find(h1); % index into data of current talker
   c1 = data(k1(1)).class; % class of current talker
   k = find(strcmp(c1,{data.class}) & ~h1); % index into data, same CLASS, different TALKERs
   
   % Pool data:
   x = [data(k(1)).features]; % scalar, ROW vector, or matrix of ROW vectors, grow dynamically
   for p1=2:length(k),
      x = [x;data(k(p1)).features]; % stack scalars or ROWs of features
   end;
   
   % Train model:
   modelsLOO(p).obj = gmdistribution.fit(x,parameters.numMixtures,'covtype',parameters.covType,...
      'sharedcov',parameters.sharedCov,'options',GMMoptions);
   modelsLOO(p).talker = talkerID{p};
end;

% Init confusion matrices for test results on all samples of a talker, ground truth in each COLUMN,
% model output in each ROW:
cmAllSample = zeros(length(modelsAll)); % test-on-train experiment
cmLOOSample = zeros(length(modelsAll)); % leave-one-out experiment

% Test samples from each talker:
for p=1:length(talkerID), % for each test talker
   % Get indices of test samples:
   h1 = strcmp(talkerID{p},{data.talker}); % logical output, case sensitive
   k1 = find(h1); % index into data of current talker
   c1 = data(k1(1)).class; % class of current talker
   c1Index = find(strcmp(c1,classID)); % index into classID, ground truth
   
   % Init confusion matrices, ground truth in each COLUMN, model output in each ROW:
   cmAll = zeros(length(modelsAll)); % test-on-train experiment
   cmLOO = zeros(length(modelsAll)); % leave-one-out experiment
   
   % Test each sample from current talker:
   LLSample = zeros(length(k1),length(modelsAll)); % log-likelihood for all samples, test-on-train experiment
   LLLOOSample = zeros(length(k1),length(modelsAll)); % log-likelihood for all samples, leave-one-out experiment
   for p1=1:length(k1),
      % Get current sample:
      x = data(k1(p1)).features;
      
      % Get likelihood from modelsAll:
      LL = zeros(1,length(modelsAll));
      for p2=1:length(modelsAll), % for each model
         % Compute likelihood of x for current model:
         xLik = pdf(modelsAll(p2).obj,x);
         
         % Compute log-likelihood:
         LL(p2) = sum(log(xLik),1);
      end;
      
      % Get likelihood from modelsLOO:
      xLikLOO = pdf(modelsLOO(p).obj,x);
      LLLOO = LL; % copy LL from all models
      LLLOO(c1Index) = sum(log(xLikLOO),1); % replace with LOO log-likelihood
      
      % Update LLSample:
      LLSample(p1,:) = LL;
      LLLOOSample(p1,:) = LLLOO;
      
      % Determine winning models:
      [junk,mAllIndex] = max(LL); % index of winning model, test-on-train experiment
      [junk,mLOOIndex] = max(LLLOO); % index of winning model, leave-one-out experiment
      
      % Update confusion matrices:
      cmAll(mAllIndex,c1Index) = cmAll(mAllIndex,c1Index)+1;
      cmLOO(mLOOIndex,c1Index) = cmLOO(mLOOIndex,c1Index)+1;
   end;
   
   % Determine winning model using all samples:
   [junk,mAllSampleIndex] = max(sum(LLSample,1)); % index of winning model, test-on-train experiment
   [junk,mLOOSampleIndex] = max(sum(LLLOOSample,1)); % index of winning model, leave-one-out experiment

   % Update confusion matrices:
   cmAllSample(mAllSampleIndex,c1Index) = cmAllSample(mAllSampleIndex,c1Index)+1;
   cmLOOSample(mLOOSampleIndex,c1Index) = cmLOOSample(mLOOSampleIndex,c1Index)+1;
   
   % Save results:
   output(p).cmAll = cmAll;
   output(p).cmLOO = cmLOO;
   output(p).accuracyAll = sum(diag(cmAll))/length(k1)*100; % percent correct
   output(p).accuracyLOO = sum(diag(cmLOO))/length(k1)*100; % percent correct
   output(p).talkerID = talkerID{p};
   output(p).classID = c1Index;
   output(p).classAllSample = mAllSampleIndex;
   output(p).classLOOSample = mLOOSampleIndex;
end;

% Tabulate confusion matrices over all talkers:
cmAll = output(1).cmAll;
cmLOO = output(1).cmLOO;
for p=2:length(output), % for each talker
   cmAll = cmAll+output(p).cmAll;
   cmLOO = cmLOO+output(p).cmLOO;
end;

return;

function [classID,talkerID,data] = getDataList(data)
% This function returns a list of classes and talkers in data. The class and talker IDs may be
% scalars or strings.  All scalars are converted to strings.

% Check that class and talker labels are strings, convert if necessary:
for p=1:length(data),
   data(p).class = convertDataEntry(data(p).class);
   data(p).talker = convertDataEntry(data(p).talker);
end;

% Get unique list:
classID = unique({data.class});
talkerID = unique({data.talker});

return;

function y = convertDataEntry(x)
% This function converts x to a string y.  x may be a cell array, numeric, or logical value.

% Check for cell array entry, store first element only if necessary:
if iscell(x),
   x = x{1};
end;

% Convert x to string:
if ischar(x),
   y = x;
else % not string, so convert
   if isnumeric(x),
      y = num2str(x);
   elseif islogical(x),
      y = num2str(double(x));
   end;
end;

return;

% Bye!
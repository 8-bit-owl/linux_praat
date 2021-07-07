function output = calculateHFCCmeasures(X)
% This function calculates various measures for matrix X (cepstral coefficients, delta coefficients,
% or delta-delta coefficients).
% Input:
%    X -- real NxL matrix, N cepstral coefficients, L frames
% Output:
%    output -- struct of output measures
%       .stDevVector -- Nx1 real vector, standard deviation of each ROW of X
%       .stDevSum -- real scalar, standard deviation along each ROW of X, summed
%       .varSum -- real scalar, variance along each ROW of X, summed
%       .euclDistMean -- real scalar, average of Euclidean distance between adjacent COLUMNS of X
%       .absDistMean -- real scalar, average of L1 distance between adjacent COLUMNS of X

% Mark Skowronski, May 29, 2013

% Init output:
output = struct([]);

% Standard deviation vector and sum:
output(1).stDevVector = std(X,[],2); % standard deviation along each ROW
output(1).stDevSum = sum(output.stDevVector); % standard deviation along each ROW

% Variance sum:
output(1).varSum = sum(output.stDevVector.^2); % variance along each ROW

% Euclidean distance average:
DX = abs(X(:,2:end)-X(:,1:end-1)); % absolute difference matrix
output(1).euclDistMean = mean(sqrt(sum(DX.^2,1)));

% Absolute distance average:
output(1).absDistMean = mean(sum(DX,1));

return;

% Bye!
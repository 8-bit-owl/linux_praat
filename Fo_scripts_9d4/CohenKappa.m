function [Kexact,Kadjacent] = CohenKappa(x,y)
% This function calculates Cohen's kappa for pairs of ratings x and y. Kappa is calculated for exact
% and adjacent (plus-minus one rating point) agreement.
% Input:
%    x -- Nx1 integer vector, N ratings in range [1,D] from rater 1
%    y -- Nx1 integer vector, N ratings in range [1,D] from rater 2, paired with ratings in x
% Output:
%    Kexact -- real scalar, kappa for exact agreement
%    Kadjacent -- real scalar, kappa for adjacent agreement
% Kappa = (#correct - #chance)/(#total - #chance)
%    #total = N
%    #correct = sum(abs(x-y)<=0) for exact, sum(abs(x-y)<=1) for adjacent agreement
%    #chance: calculated from confusion matrix tabulated from x and y

% Tabulate confusion matrix between x and y:
N = length(x);
D = max(max(x),max(y)); % largest rating in either input
cm = zeros(D);
for p=1:N,
   i = x(p);
   j = y(p);
   cm(i,j) = cm(i,j)+1;
end;

% Tabulate chance matrix from cm:
s1 = sum(cm); % sum of each COLUMN
s2 = sum(cm'); % sum of each ROW
ch = (s2'*s1)/N; % ch(i,j) is chance number of times rater 1 responded "i" and rater 2 responded "j".

% Calculate exact and adjacent counts:
numCorrectExact = sum(diag(cm,0));
numChanceExact = sum(diag(ch,0));
numCorrectAdjacent = sum(diag(cm,-1))+sum(diag(cm,0))+sum(diag(cm,1));
numChanceAdjacent = sum(diag(ch,-1))+sum(diag(ch,0))+sum(diag(ch,1));

% Calculate exact and adjacent kappa:
Kexact = (numCorrectExact - numChanceExact)/(N - numChanceExact);
Kadjacent = (numCorrectAdjacent - numChanceAdjacent)/(N - numChanceAdjacent);

return;

% Bye!
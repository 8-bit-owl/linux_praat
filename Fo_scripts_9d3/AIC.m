function a = AIC(MSE,N,K)
% This function calculates the Akaike information criterion from modeling error.
% Input:
%    MSE -- real scalar, mean-squared error from a model
%    N -- integer scalar, number of samples used to train model
%    K -- integer scalar, number of parameters in model
% Output:
%    a -- real scalar, AIC value

% Mark Skowronski, March 13, 2013

a = N*log(MSE) + 2*K;

return;

% Bye!
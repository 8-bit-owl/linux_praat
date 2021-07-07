function a = AICc(MSE,N,K)
% This function calculates the corrected Akaike information criterion from modeling error, which is
% prefered over AIC when N/K < 40 (small dataset regime).
% Input:
%    MSE -- real scalar, mean-squared error from a model
%    N -- integer scalar, number of samples used to train model
%    K -- integer scalar, number of parameters in model
% Output:
%    a -- real scalar, AICc value

% Mark Skowronski, March 13, 2013

a = N*log(MSE) + 2*K*N/(N-K-1);

return;

% Bye!
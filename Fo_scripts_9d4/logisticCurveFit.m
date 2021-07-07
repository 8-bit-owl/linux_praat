function output = logisticCurveFit(x,y,IC)
% This function returns curve fit info for fitting x to y using a logistic function.
% Logistic equation:
%  y(x) = 100/(1+exp(-s*(x-x0)/25))
%     x0 -- logistic midpoint, y(x0) = 50
%     s -- midpoint slope, dy(x0)/dx = s
%  Note: The lower and upper asymptotes are assumed to be 0 and 100, respectively.
% Input:
%    x -- Nx1 real vector, independent variable
%    y -- Nx1 real vector, dependent variable
%    IC -- 1x2 real vector, [x0,s] initial conditions
% Output:
%    output -- struct of output variables
%       .fitResults -- cfit object, from fit.m
%       .gof -- struct, goodness-of-fit values for curve fit, from fit.m

% Mark Skowronski and Lisa Kopf, March 17, 2014

% Ensure x, y are COLUMN data:
x = x(:);
y = y(:);

% Define logistic equation and variables:
eq = '100/(1+exp(-s*(x-x0)/25))'; % logistic function 
dep = {'x0','s'};
ind = {'x'};

% Create model and options:
model = fittype(eq,'coefficients',dep,'independent',ind);
fitOptions = fitoptions('method','nonlinearleastsquares','startpoint',IC);

% Fit model to data:
[output(1).fitResults,output(1).gof] = fit(x,y,model,fitOptions);
output(1).fitResults

return;

% Bye!
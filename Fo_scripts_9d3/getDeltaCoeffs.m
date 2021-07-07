function Dx = getDeltaCoeffs(x,deltaSize)
% This function calculates delta coefficients across each ROW for the matrix in x.  +/- deltaSize
% frames are used to calculate delta coeffs.

% Construct delta filter:
h = [deltaSize:-1:-deltaSize];

% Filter each ROW in x:
Dx = conv2([1],h,x,'same'); % convolve across rows w/ [1], then across columns with h, return same size as x

return;


function y = etaprime(x, threshold)
% ETAPRIME is the derivative of the soft threshold function. In reality it
% returns 0 if x in [-threshold, threshold] and 1 otherwise.

y = (x > threshold) + (x < -threshold);
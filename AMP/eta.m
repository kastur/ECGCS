function y = eta(x, threshold)
% ETA performs a soft thresholding on the input x.

y = wthresh(x, 's', threshold);
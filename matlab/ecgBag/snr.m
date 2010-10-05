function ratio = snr(sig, noise);
% ratio = snr(sig, noise);
%
% GC 2005

ratio = 10*log10(var(sig)/var(noise));

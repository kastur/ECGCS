function [Px, Prob] = lomb(t, x, f)
%   [Px Prob] = lomb(t, x, freq)
%
%   Calculates the Lomb-Scargle normalized periodogram values 
%   "Px" as a function of the supplied vector of frequencies 
%   "freq" for input vectors "t" (time) and "x" (observations).
%   Also returns the probability "Prob" that the null hypothesis
%   is valid (same length as Px and freq). Time stamps, t and 
%   amplitudes "x" must be the same length. 
%
%   See Scargle J.D.:"Studies in astronomical time series analysis. II. 
%   Statistical aspects of spectral analysis of unevenly spaced data,"
%   Astrophysical Journal, vol 263, pp. 835-853, 1982.  ... and
%   Lomb N.R: "Least-squares frequency analysis of unequally spaced data",
%   Astrophysical and Spcae Science, vol 39, pp. 447-462, 1976.

% This file is an adaptation of the mysterious lomb.m which was emailed
% to me some time ago by a colleague who obtained it from a forgotten source.
% I claim no responsibility for its accuracy, although it seems to
% correspond with lomb.c from NRinC (but not fasper.c)
% Any information you have on this file, please email me:
% gari AT physionet DOT org


if nargin < 2
 error('must have an amplitude for each time stamp')
end

% If no frequency vector is supplied, invent a default up to the 
% highest frequency available (> Average Nyquist)
if nargin < 3
 maxfreq = 1/min(diff(t));
 f = [1/512:1/512:maxfreq];
end

% check length of inputs
if length(t) ~= length(x); 
 error('t and x not same length');	
 exit; 
end;

% subtract mean, compute variance, initialize Px
z   = x - mean(x);
var = std(x);
N   = length(f);
Px  = zeros(size(f));

% compute power by looping over all frequencies
for i=1:length(f)
    w=2*pi*f(i);
    if w > 0 
       twt = 2*w*t;
       tau = atan2(sum(sin(twt)),sum(cos(twt)))/2/w;
       wtmt = w*(t - tau);
       Px(i) = (sum(z.*cos(wtmt)).^2)/sum(cos(wtmt).^2) + ...
		(sum(z.*sin(wtmt)).^2)/sum(sin(wtmt).^2);
     else
	Px(i) = (sum(z.*t).^2)/sum(t.^2);
     end
end


% normalize by variance and compute probabilities at each frequency

for i=1:length(Px)          
    if var~=0            % check for divide by zero
        Px(i)=Px(i)/2/var.^2;
        Prob(i) = 1-(1-exp(-Px(i)))^N;
    else
        Px(i)=inf; 
        Prob(i)=1;
    end;
    if Prob(i) < .001  % allow for possible roundoff error
	Prob(i) = N*exp(-Px(i));
    end
end;


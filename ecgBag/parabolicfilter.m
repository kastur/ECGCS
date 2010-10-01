function [baseline, corrected, coeffs] = parabolicfilter(x, n)

% [baseline corrected coeffs] = parabolicfilter(x, n); 
% designs an n-point parabolic filter according to 
% c = -3*(7 + 20*j.*j - 3*n^2))/(4*(n^3 -4*n)); where j runs from 
% -delay to + delay (delay = (n-1)/2 and must be odd). Zero phase 
% filtering is then performed with filtfilt.m to produce a baseline:
% baseline = filtfilt(coeffs, 1, y);
% corrected = x-baseline and coeffs are the filter coefficients used.
% In other words, for each point in x, the filter fits a parabola 
% to the (n-1)/2 preceding points, the current point, and the (n-1)/2 
% following points, and is evaluated at the current point.
% Based on some old medical code converted by David Clifton into C,
% rewritten by Shelley Cazares into Matlab in 1999 and updated by
% Gari Clifford in 2003. 

% assume a sampling frequency of 256Hz and we are filtering BP signal
if nargin<2
    n = round(256*1.3);
end
% check the window size is odd
if(rem(n,2)==0)
    error('n must be odd')
end
%%%%% Length of x.                                        
L = length(x);        

%%%%% Delay of filter.
delay = (n-1)/2;      

%%%%% Coefficients of filter.
coeffs = zeros(n+1,1);       
coeffsindex = [-delay:1:delay]';
coeffs = (-3*(7 + 20*coeffsindex.*coeffsindex - 3*n^2))/(4*(n^3 -4*n));

%%%%% Pad x with (n-2)/2 zeros at the end.  Matlab with automatically
%%%%% pad x with (n-1) zeros at the beginning.  Filter works on chunks
%%%%% of data including the current sample and the (n-1) preceding samples.
%%%%% Completely causal.
padx = zeros(L+delay,1); 
padx(1:L) = x;           

%%%%% Filter x.
shifted = zeros(L+delay,1);              % Initialize the filtered data.
shifted = filter(coeffs, 1, padx);       % Filter the data.
y = shifted(delay+1:delay+L);            % Compensate for the (n-1)/2 delay.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Due to the (n-1)/2 delay, the first and last (n-1)/2 filtered points are
% not a good fit to the data.  This is because the first and last (n-1)/2
% points are fit to at least one of the padded zeros at the beginning and
% end of the data.  Therefore, for the first and last (n-1)/2 points, we will
% output the unfiltered data, instead of the filtered data.  For all other
% points, we will output the filtered data.

y(1:1+delay) = x(1:1+delay);
y(L-delay:L) = x(L-delay:L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseline = filtfilt(coeffs, 1, y);
corrected = x-baseline;

function [real_t,real_R] =  down_sample_realistically(x,t,samp_freq,noise_level)

% [real_t,real_R] =  down_sample_realistically(x,t,samp_freq,noise_level);
% noise level is a % - so 100/length(x) is 1 point!
% Defaults are: samp_freq=7Hz, noise level=0;

% These ECG processing routines are made available under the 
% GNU general public license. If you have not received a copy 
% of this license, please download from http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% (C) G. Clifford 2001 : gari AT ieee DOT org 



format long e

% check if defaults are changed
if nargin < 4
   noise_level=0
end

if nargin < 3
   samp_freq=7
end

if nargin < 2
   error('You need at least a time vector and a signal')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_points = length(x);
i=1;
q=1;
real_R(1)=x(i);
real_t(1)=t(i);
last_t=t(i);
for i=2:no_points;
  if( (t(i)-last_t) > x(i) )  % then it's good
     q=q+1;     
     real_R(q)=x(i);
     real_t(q)=t(i);
     last_t   =t(i);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(t,x,'.b'); hold on; plot(real_t,real_R,'+--r')


function [rr, rs] = rrextract(ecg, sf, thresh, testmode, twind)
% [rr rs] = rrextract(ecg, sf, thresh, testmode, twind) extracts rr times 
% and rs amplitudes from ecg with sampling frequency sf using rpeakdetect.m 
% (see RPEAKDETECT for a detailed explanation of how to set the threshold
% value 'thresh'). This wrapper function calls the batch process 
% RPEAKDETECT on every twind seconds (default 1000) and strings the results 
% together. If testmode is set to ==1 then rpeak detect will pause every
% twind seconds and plot 4 graphs - 
% 1) ECG with zero phase band passed (BP) filtered ECG over the top of it
% 2) ECG with differentiated, squared & integrated (energy) ECG overlayed 
%    on BP filtered ECG (together with green and magenta circles where 
%    thresh*max(energyECG) is breached; defining the peak search region)
% 3) Band passed filtered ECG with R & S points marked (max and min of 
%    BP filtered ECG inside search region)
% 4) RR interval time series
% Currently only BP filters for 128 and 256Hz are implemented. If you feed
% in a signal with any other sf, no BP filtering is done. If you want 
% to add filter support for other sf's then the load in the coefficients in 
% filterECG256Hz.m into sptool, change the frequency, and output the new
% coeffs. Please try to follow the format of filterECG256Hz - FIR, odd no.
% of taps, check the filter perfomance then send me the coeffs and I'll
% add support for the new sf. Note that twind must be more than 3 times 
% as long as the longest filter order for filtfilt.m to work (currently 
% 313 taps for the 256Hz filter and 245 taps for the 128Hz filter).
% e.g. [rr rs] = rrextract(data(:,1), 256, 0.2, 0, 1000); 
% is the same as [rr rs] = rrextract(data);
%
% rpeakdetect.m was written by G.D.Clifford (gari@ieee.org) and the 
% rrextract.m wrapper was written by P.E.McSharry and G.D.Clifford. 
% These ECG processing routines are made available under the 
% GNU general public license. If you have not received a copy 
% of this license, please download from http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% gari@ieee.org

if nargin < 5
   twind = 1000;
end
%%%%%%%%%%% testmode == 1 indicates we should graph the data
if nargin < 4
   testmode = 0;
end
%%%%%%%%%%% make threshold default 0.2 -> this was 0.15 on MIT data 
if nargin < 3
   thresh = 0.2;
end
%%%%%%%%%%% make sample frequency default 256 Hz 
if nargin < 2
   sf = 256;
end
%%%%%%%%%%% make sample frequency default 256 Hz 

if (twind<11)
     error('window length must be at least 11 seconds')
end

% check dimension of ecg
if(min(size(ecg))>1)
     error('ecg should be 1 dimensional');
end

% set length of window
W = twind*sf;
N = length(ecg);

Nw = floor(N/W);


i2 = 0;
for i=1:Nw
   data = ecg((i-1)*W+1:i*W,:);
   [hrv R_t R_amp R_index S_t S_amp]  = rpeakdetect(data,sf,thresh,testmode);
   Nr = length(R_index);
   i1 = i2;
   i2 = i2+Nr; 
   rr(i1+1:i2) = (i-1)*W + R_index;
   rs(i1+1:i2) = R_amp-S_amp;
end
data = ecg(Nw*W+1:end,:);
[hrv R_t R_amp R_index S_t S_amp]  = rpeakdetect(data,sf,thresh,testmode);
Nr = length(R_index);
i1 = i2;
i2 = i1+Nr;
rr(i1+1:i2) = Nw*W + R_index;
rs(i1+1:i2) = R_amp-S_amp;

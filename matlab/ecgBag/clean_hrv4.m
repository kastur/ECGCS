function [hrv] = clean_hrv4(temp_hrv, pc)

% [hrv] = clean_hrv4(temp_hrv, pc)
% Takes the two col (t, RR) HRV data (temp_hrv) and removes any 
% RR intervals <> 1-pc% of the previous (`true' recorded) interval. 
% Default value is 80% (== removal of beats <> 20% change over previous)
% See Clifford G.D., McSharry P.E., Tarassenko L.: "Characterizing 
% Artefact in the Normal Human 24-Hour RR Time Series to Aid Identification 
% and Artificial Replication of Circadian Variations in Human Beat to Beat
% Heart Rate Using a Simple Threshold", Computers in Cardiology, 2002.
%
% These ECG processing routines are made available under the 
% GNU general public license. If you have not received a copy 
% of this license, please download from http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% G. Clifford : gari@ieee.org


if nargin < 2
   pc=80;
end

if nargin < 1
   error('must have at least one argument')
end


%% Initialise
j=1;

%% check validity of first sample
%% (assuming most data is good and fairly stationary.)
mean_tmp=mean(temp_hrv(:,2)); 
stp=0;             %
st=size(temp_hrv); %
if st(1)>1         % number of rows must be > 1
    while ((stp==0) & ((abs(temp_hrv(j,2)-mean_tmp)>10 )&(abs(temp_hrv(j,2)-temp_hrv(j+1,2)) > 10)))
        j=j+1;
        if j==st(1)% if no row is found that fits criteria,
            stp=1; % exit loop
        end;       %
    end;
end;    
%% if the first j samples are not good, skip them
%% make first valid hrv value.
hrv=[];
if stp==0    % only if good samples had been found
i=1;
hrv(i,1)=temp_hrv(j,1);
hrv(i,2)=temp_hrv(j,2);
valid_RR = temp_hrv(j,2);
valid_T  = temp_hrv(j,1);
sz=length(temp_hrv);

% a bodge to allow the loop to run over the end 
%(mirror around the last sample by one sample)
temp_hrv(sz+1,2)=temp_hrv(sz-1,2);
% and add the correct time stamp - last time stamp plus new RR value
temp_hrv(sz+1,1)=temp_hrv(sz,2)+temp_hrv(sz+1,2);

% set hrv counter.
k=2;
l=0;
for i=1:sz-j; % count through all data points starting from 
              % one after the first valid one
              % and keep only if it doesn't change by more than pc percent
              % on the last beat or on the next beat                 

  if (abs((60/temp_hrv(j+i,2))-(60/hrv(k-1,2))) < ((1-( (pc-(l*3))/100))*(60/hrv(k-1,2)))) | (abs((60/temp_hrv(j+i,2))-(60/temp_hrv(j+i+1,2))) < ((1-( (pc-(l*3))/100))*(60/hrv(k-1,2)))) 
       hrv(k,1)=temp_hrv(j+i,1);
       hrv(k,2)=temp_hrv(j+i,2);
       k=k+1; % counter for each beat kept
       l=0; % reset local counter if we keep a beat
      else
       l=l+1; % local counter for consecutive beats removed (for each beat removed we decrement pc by (l*3)% ... see (pc-(l*3)) in if statement above
      end
end
end;

%save cleaned_hrv2.dat hrv -ascii

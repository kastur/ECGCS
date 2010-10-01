function [dataout] =  select_train(ecgdata, fiddata, N_train, n_lead, downsample, wind, offset, samp_freq)

% [dataout] =  select_train(ecgdata, fiddata, N_train, n_lead, downsample, wind, offset, samp_freq)
% pass in an n_lead ecg and form a set of training vectors 
% from this ecg data
% ecgdata - raw ecg - timing
% fiddata - indexes of the fiducial points
% N_train - number of training vectors
% n_lead - the lead you want
% samp_freq - sampling frequency
% downsample - training vector downsample factor
% wind - window over which you want to segment training vectors
% offset - number of samples offset at start of ECG

% m = size of window needed around QRS ... in_sz = round(samp_freq/m)*2 +1;
% ... becoming redundant soon!

m = 3;

% check if defaults are changed
if nargin < 8
   samp_freq = 256; %smp
end

if nargin < 7
   offset = ecgdata(1,1);
end

if nargin < 6
   wind = 0.5; % seconds
end

if nargin < 5
   downsample = 4;
end

if nargin < 4
   n_lead = 1; %smp =2  - - MIT ==1
end

if nargin < 3
   N_train =  length(fiddata);
end

if nargin < 2
   error('Must have data and fiducial points.');
end


%offset_ind = round( (offset+5*(1/samp_freq))*samp_freq ); 
offset_ind = round( offset*samp_freq ); 

in_sz = round(samp_freq/m)*2 +1; % no of input nodes

%dataout = zeros(N_train,in_sz); % intialise the training matrix 

for i=1:N_train % collect N_train QRS complexes
                                      % approximate the sample point
  % actually - don't - just do it directly
  fid_index=fiddata(i);  %round(samp_freq*fiddata(i))-offset_ind; 
  %size(dataout) - debug
  if (downsample == 0) % don't downsample - do it quick
    temp = ecgdata(fid_index-round(samp_freq/m):fid_index+round(samp_freq/m),n_lead); 
    dataout(i,1:length(temp)) =  temp'; % put in training vector to output matrix
   
  else 
    len_array=round((wind*samp_freq)/(downsample));


% fid_index % debug
    % put in points below fid point
    for j=1:round(len_array/2)  
    %  round(len_array/2)-j+1
      temp(round(len_array/2)-j+1) = ecgdata(fid_index-(j*downsample),n_lead);
    end
    temp(j+1) = ecgdata(fid_index,n_lead);  % middle point
    %j+1
    % put in points above fid point
    for k=1:round(len_array/2)
      temp(j+k+1) = ecgdata(fid_index+(k*downsample),n_lead);
      %j+k+1
    end
   dataout(i,1:length(temp)) =  temp; % put in training vector to output matrix 
  
  end

end



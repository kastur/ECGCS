function [dataout] =  recon_ECG(ecgorig, ecgmesh, fiddata, N_train, wind, samp_freq)

% [dataout] =  recon_ECG(ecgorig, ecgmesh, fiddata, N_train, wind, samp_freq);
% pass in a (filtered) mesh ECG and an original ECG and fiducial point data. 
% unwrap to the filtered mesh of QRS complexes and add to the original 
% ecg at each fiducial point.
% ecgorig - raw ecg 
% ecgmesh - filtered mesh of QRS complexes
% fiddata - indexes of the fiducial points
% N_train - number of training vectors
% wind - window width of each QRS complex
% samp_freq - sampling frequency


% m = size of window needed around QRS ... in_sz = round(samp_freq/m)*2 +1;
% ... becoming redundant soon!

m = 3;

% check if defaults are changed
if nargin < 6
   samp_freq = 256; %smp
end
if nargin < 5
   wind = 0.5; % seconds
end
if nargin < 4
   N_train =  length(fiddata);
end
if nargin < 3
   error('Must have original data, filtered mesh data and fiducial points.');
end

downsample=1;

%offset_ind = round( (offset+5*(1/samp_freq))*samp_freq ); 
% offset_ind = round( offset*samp_freq ); 

in_sz = round(samp_freq/m)*2 +1; % no of input nodes

% intialise data out with original data
dataout=ecgorig;

len_array=round((wind*samp_freq)/(downsample));
halfwin=floor(len_array/2);

for i=1:N_train % loop over the number of QRS complexes we are replacing
f_ind=fiddata(i);
a=f_ind-(halfwin);
b=f_ind+halfwin;
%length(ecgmesh(i,:))
%fprintf('%i : %i\n',a,b);
dataout(a:b)=ecgmesh(i,:);
end

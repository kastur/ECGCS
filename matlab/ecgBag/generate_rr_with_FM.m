function [t,x] =  generate_rr_with_FM(LF,HF,HR,len,LF_amp,HF_amp,samp_freq,sigma,lf_amp_range,hf_amp_range,phi);

% [t x] =  generate_rr_with_FM(LF,HF,HR,len,LF_amp,HF_amp,samp_freq,sigma,lf_amp_range,hf_amp_range,phi);
% returns a time vector and column of evenly (samp_freq) sampled 
% RR intervals  for 'len' seconds generated from two sines at LF and 
% HF Hz. Their amplitudes are modulate by a random amount about the LF_amp
% and HF_amp with ranges 'lf_amp_range' and 'hf_amp_range'. Their frequencies
% slowly increase between LF_lo = 0.04 and LF_hi = 0.15 AND HF_lo = 0.15 
% and HF_hi = 0.4, dwelling in each frequency bin for an amount of
% time dictated by a normal distribution with SD sigma.
% phi is the fraction of 2*pi that the phase will drift between LF and HF
% over len seconds. HR is the baseline offset in BPM. 
% ** x is returned as an RR tachogram. **
% note that the quantisation of the time in each freq bin means t>x 
%
% Defaults LF = 0.095Hz, HF = 0.275Hz, HR = 60BPM, phi= 0.,
% samp_freq=7Hz, len=300s, LF_amp=+/-2BPM, HF_amp=+/-2.5BPM, 
% sigma=10000; lf_amp_range=0; hf_amp_range=0;
%
% Written by G. Clifford (C) Oxford University 2001

format long e

% check if defaults are changed
if nargin < 11  
   phi=0 % no drift
end

if nargin < 10
   hf_amp_range=0 % no AM
end

if nargin < 9
   lf_amp_range=0 % no AM
end

if nargin < 8
   sigma=1 %10000 % was 3 .. but 10000 gives no FM
end

if nargin < 7
   samp_freq=7
end

if nargin < 6
   HF_amp=2.5
end

if nargin < 5
   LF_amp=2
end

if nargin < 4
   len=300
end

if nargin < 3
   HR=60
end

if nargin < 2
   HF=0.275
end

if nargin < 1
   LF=0.095
end

%%%%%%%%%% initialise

%lf_amp_range= (1/2);    % +/- 0.5  beat per min == 1BPM
%hf_amp_range= (1.25/2); % +/- 0.625 beat per min == 1.25BPM

% how many different frequencies should we have?
% one for each sample point!
Num_of_W=round(len*samp_freq);

LF_lo = 0.04;
LF_hi = 0.15;
HF_lo = 0.15; 
HF_hi = 0.4;
centre_freq_hi=0.275;
centre_freq_lo=0.095;
sigma_hi=(HF_hi-HF_lo)/sigma;
sigma_lo=(LF_hi-LF_lo)/sigma;

fh_inc = 0.000125;
Num_of_samples=round(len*samp_freq);
hf_div= (HF_hi-HF_lo)/Num_of_W;
lf_div= (LF_hi-LF_lo)/Num_of_W;
wh_inc = 2*pi*hf_div; 
wl_inc = 2*pi*lf_div; 
wl=0;
wh=0;

% create the time vector
t= ([1:round(len*samp_freq)])/samp_freq;

% calculate the number of frequencies in each band

for (j=1:Num_of_W)
  Wh(j)= (HF_lo + ((j)*hf_div))*2*pi;
  no_samps_h(j)=Num_of_samples*(1/sqrt(2*pi*sigma^2))*exp((-0.5)*(((Wh(j)-(2*pi*centre_freq_hi))^2)/((sigma_hi)^2)));
end

normalisation = sum(no_samps_h)/Num_of_samples;
no_samps_h=no_samps_h/normalisation;
no_samps_h=round(no_samps_h);

figure;
stem(Wh/(2*pi),no_samps_h,'.')

if(LF_amp>0)
 for (j=1:Num_of_W)
   Wl(j)= (LF_lo + ((j)*lf_div))*2*pi;
   no_samps_l(j)=Num_of_samples*(1/sqrt(2*pi*sigma^2))*exp((-0.5)*(((Wl(j)-(2*pi*centre_freq_lo))^2)/((sigma_lo)^2)));
 end
normalisation = sum(no_samps_l)/Num_of_samples;
no_samps_l=no_samps_l/normalisation;
no_samps_l=round(no_samps_l);
hold on
stem(Wl/(2*pi),no_samps_l,'r.')
end

% fraction of 2*pi that phase will drift by over len seconds = phi
% phase will change by ??? radians each sample
phi_diff = (2*pi*phi)/(len*samp_freq);

% if there is no incremental frequency change - 
wh = 2*pi*centre_freq_hi;
wl = 2*pi*centre_freq_lo;

% if there is no incremental amplitude change - 
lf_amp_add=0;
hf_amp_add=0;

k=1;
for(i=1:Num_of_W) % for each decimation of the frequency
    if (no_samps_h(i)>0) % ignore zero entry indices
      for(j=1:no_samps_h(i)) % how many samples for this frequency?
	  phi_diff*(k-1);
	  lf_amp_add=lf_amp_range*2*(rand-0.5);
	  hf_amp_add=hf_amp_range*2*(rand-0.5);	  
	  wh  = (2*pi*HF_lo)+(i*wh_inc); %  
	  wl  = (2*pi*LF_lo)+(i*wl_inc); %  
	  X(k)=HR+(HF_amp+hf_amp_add)*sin(wh*t(k))+(LF_amp+lf_amp_add)*sin((wl+(phi_diff*(k-1)))*t(k));
	  k=k+1; 
      end
    end
end               % for each decimation of the frequency
 

%artificial RR - convert to RR intervals
x = 60./X;
t = [1:1:length(x)]/samp_freq;

xlabel('Frequency (Hz)')
ylabel('Number of samples at a given frequency')
hold off

% t and x are going to be unequal lengths!!!
% [Px,Fx] = psd(xx(1:2048),2048,samp_freq,hamming(2048),'mean');
% plot(Fx,20*log10(Px))       
% axis([0 0.5 -110 0])
% xlabel('Frequency (Hz)')
% ylabel('Power (dB)')
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



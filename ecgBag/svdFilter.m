function [svdfilt, resid, noise]  = svdFilter(ecg,nEVs,Fs,wind,n,thresh,noiseAmp)

% [svdfilt resid noise]  = svdFilter(ecg,nEVs,Fs,wind,n,thresh,noiseAmp);
% 
% Filter a single channel of ECG by peak detecting,
% using a threshold of thresh [default 0.2],
% aligning beats at R-peaks, performing a truncated
% SVD on the matrix of beats with n components, 
% [default n=1], reconstructing the truncated 
% decomposition matrix, then unwrapping the matrix and
% placing back into the unfiltered ECG at each of the 
% on the time time using the correct R-peak fiducual points.
% into the vector svdfilt with a residual difference of resid.
% normally distributed noise is added with noiseAmp [default=0]
% up to same variance as ecg (noiseAmp=1); 
%
% nEVs= number of eigenvalues
% Fs = sampling frequency of input signal
% wind = window in seconds to segment around ECG 
% n = 
% thresh = QRS detector threshold = 0.2
% noiseAmp = noise amplitude if additive noise if required.
%
% GNU public license, (C) Gari D Clifford 2004 gari AT physionet DOT org


if nargin<7
noiseAmp=0;
end
if nargin<6
thresh=0.2;
end
if nargin<5
n=1;
end
if nargin<4
wind=1000;    
end
if nargin<3
Fs=256;
fprintf('Assuming sampling frequency of %f Hz\n',Fs);
end
if nargin<2
nEVs=5;
end

noise=0;

if exist('ecg')==0 | length(ecg)<1000
fprintf('ecg not defined or too small, generating one myself thankyou!');
[ecg, ipeaks] = ecgsyn24(Fs,Fs,0,0,0,0,256,60,1,0.5,0);
%[ecg, ipeaks] = ecgsyn24(Fs,Fs,rr,ti,ai,bi,N,hrmean,hrstd,lfhfratio,Anoise,seed);
%[rr] = rrprocess(Fs,Fs,0,60,1,0.5,256);

% add noise to ecg?
noise = randn(size(ecg));
noise=noise/std(noise); % unit variance
noise = noise*std(ecg);
noise-noise*noiseAmp;
ecg = ecg+noise;
end


j=min(size(ecg));
[rr rs] = rrextract(ecg(:,j), Fs, 0.2, 0, 1000); 

if(wind>100)
 wind=mean(diff(rr))/Fs;
end

[meshecg] = select_train(ecg(:,j), rr(2:length(rr)-1), length(rr)-2, 1, 1, wind, 0, Fs);

[U,S,V] = svds(meshecg,nEVs);

recon = U*S*V'; %' the decomposition of the data IS .... U S V^T 

[svdfilt] =  recon_ECG(ecg, recon,rr(2:length(rr)-1), length(rr)-2, wind,Fs);

resid=ecg-svdfilt; 

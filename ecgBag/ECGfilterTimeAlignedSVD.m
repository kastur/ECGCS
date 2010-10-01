function filtecg = ECGfilterTimeAlignedSVD(ecg, Fs, wind, QRSthresh, evS)
% filtecg = ECGfilterTimeAlignedSVD(ecg, Fs, wind, QRSthresh, evS)
% 
% use QRS detector and SVD to filter ECG
% 
% e.g.
% filtecg = ECGfilterTimeAlignedSVD(ecg, 256, 1, 0.3);
% 
% G D Clifford 2004 gari AT mit DOT edu

% SVD analysis
addpath /home/gari/CODE/TOOLS4ECGANALYSIS/
%load X
%ecg = X(1,:);

if nargin <2
 Fs=256;
end

if nargin < 3
 wind=1; % how many seconds/2 to segment either side of fiducial point 
end

if nargin < 4
 QRSthresh = 0.1;
end

if nargin < 5
 evS = 1;
end

% QRS detect
[rr rs] = rrextract(ecg', Fs, QRSthresh); %'
figure; 
plot(ecg); 
hold on;
plot(rr,ecg(rr),'+m')
title('data with fiducual points')

% make a mesh 
% [dataout] =  select_train(ecgdata, fiddata, N_train, n_lead, downsample, wind, offset, samp_freq)
[meshecg] = select_train(ecg', rr(2:length(rr)), length(rr)-1, 1, 1, wind, 0, Fs); %'
figure;
mesh(meshecg);
xlabel(' sample number')
ylabel(' beat number')
zlabel('  Amplitude')
ttl=strcat('Data realigned by fiducual points in centre of window, ',num2str(wind),'s wide')
title(ttl);

% SVD filter
[U,S,V] = svd(meshecg);
stem(diag(S))
%axis([0.5 8.5 0 12])
title('Eigenspectrum')

[U,S,V] = svds(meshecg,evS);
%filtecgmesh=U*U'*meshecg; %'
filtecgmesh=U*S*V'; %'
figure;mesh(filtecgmesh)
			 title('SVD filtered data')

% reconstruct time series
figure;

[filtecg] =  recon_ECG(ecg', filtecgmesh,  rr(2:length(rr)-1), length(rr)-2, 1, Fs); %'
subplot(3,1,1);plot(ecg);ylabel('mixed ECG');
subplot(3,1,2);plot(filtecg);ylabel('SVD filtered ECG');
subplot(3,1,3);plot(ecg'-filtecg);ylabel('Residual');%'

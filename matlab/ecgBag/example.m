% this file demonstrates how to use the routines provided
% in ecgBag.tar to process a single stream of ecg called ecg128
% For interest, the signal is only 128Hz, and so we upsample
% it with an inverse FFT to 3*128Hz. The peak detection routines
% will work at any reasonable sample rate, although above about 
% 500Hz it starts to get a bit weird. Check your output - if the 
% peaks are slightly off - use maxjiggle.m to recentre them.

% Note that the details of the QRS detector are explained in 
% J. Pan \& W. Tompkins - A real-time QRS detection algorithm 
% IEEE Transactions on Biomedical Engineering, vol. BME-32 NO. 3. 1985.
% P. Hamilton \& W. Tompkins. Quantitative Investigation of QRS 
% Detection  Rules Using the MIT/BIH Arrythmia Database. 
% IEEE Transactions on Biomedical Engineering, vol. BME-33, NO. 12.1986.
%
% ... so if you choose to play with the thresholds, you'll need to 
% play with rpeakdetect.m directly. rpeakdetect.m also gives you 
% values for R-S amplitudes, and R-S times. 
%
% A warning - this QRS detector is a bit sensitive to morphology 
% and heart rate, so be careful. Always check your output - don't
% blindly trust it.
%
% This software is offered freely and without warranty under the 
% GNU public license. Please don't blame me if anything goes wrong :)
%
% gari@ieee.org

% testmode = 0 indicates no graphing
testmode = 1;

% load example file
load ecg128.mat

% give it a time axis
samplerate0 = 128;
t1 = [1/samplerate0:1/samplerate0:round(length(ecg128)/(samplerate0))];

% filter ECG
[ecg128filt]  = filterECG128Hz(ecg128);

%%
% plot data
if (testmode==2)
figure(2);hold off;
plot(t1,ecg128-mean(ecg128));
hold on
plot(t1,ecg128filt,'r');
xlabel('time (s)')
ylabel('Amplitude (8bit)')
title('600 seconds of 128Hz 8 bit ECG filtered (red) and unfiltered (blue)')
end

%%
% resample the signal in the Frequency domain at n times the original sampling frequency 

% upsample rate
n = 2; % == 256Hz (3*128)

% time the process
tstart=cputime;
ecgInterp256 = interpft(ecg128-mean(ecg128),length(ecg128)*n);
ttaken1 = cputime-tstart;
fprintf('Time taken to upsample 600 seconds of ECG from 128 to 256Hz = %f s\n',ttaken1);

% give the new ECG a time vector
t2 = [1/(n*samplerate0):1/(n*samplerate0):round(length(ecgInterp256)/(samplerate0*n))];

ecgInterp256filt = interpft(ecg128filt,length(ecg128filt)*n);
[ecg256filtInterp]  = filterECG256Hz(ecgInterp256);
figure(3); hold off
plot(t1,ecg128filt,'c');
fprintf('First we plot the original ECG before FFT upsampling \n press any key to continue\n');
hold on;
pause
plot(t2,ecgInterp256,'g');
fprintf('Then we plot the original ECG after FFT upsampling \n press any key to continue\n');
pause
plot(t2,ecgInterp256filt,'r');
fprintf('Then we plot the FFT interpolated ECG with pre filtering (at 128Hz)\n press any key to continue\n');
pause
plot(t2,ecg256filtInterp,'b');
fprintf('Then we plot the FFT interpolated ECG with post filtering\n press any key to continue\n');
pause;

fprintf('\n So it looks like filtering THEN upsampling (red line) is the correct way to do it\n');
fprintf('\n Now we will do some peak detection ... \n');


% give them a time line
tmp = cputime;
[rrs256 RS] = rrextract(ecgInterp256filt',samplerate0*n,0.2,testmode);
ttaken2 =  cputime- tmp;
if(testmode==0)
    fprintf('Time taken to calc RR intervals = %f',ttaken2);
end

% check plot
if testmode == 2
 figure(2)
 hold off;plot(t1,ecg128-mean(ecg128));hold on;plot(t2,ecgInterp256filt,'r');
 % plot fiducial points
 plot(t2(rrs256),ecgInterp256filt(rrs256),'+g');
end

% calculate RR intervals
rr256 = diff(rrs256)/(samplerate0*n);
rrt256= rrs256(2:length(rrs256))/(samplerate0*n);
rs256= RS(2:length(rrs256));

if testmode == 1
 figure(1);
 subplot(4,1,4)
 hold off;
 plot(rrt256,rr256,'+r');
 hold on
end

% clean data
temp_hrv256 = [rrt256', rr256', rs256'];

% remove the ectopic beats and artefacts
[hrv256] = clean_hrv4([rrt256', rr256'], 88);
fprintf('%i RR interval(s) removed from timing\n',length(temp_hrv256)-length(hrv256));


[hrv256] = clean_RR_RS(temp_hrv256, 88,75);
fprintf('%i RR interval(s) removed from timing and morphology\n',length(temp_hrv256)-length(hrv256));

figure(1);
subplot(4,1,4);hold on;
plot(hrv256(:,1),hrv256(:,2),'+b');

% check PDSs of each
freq_vect =  [1/1024 : 1/1024 : 512 * 1/1024];
[Px, Prob] = lomb(hrv256(:,1), hrv256(:,2)-mean(hrv256(:,2)), freq_vect);

if testmode == 1
 figure(2);
 hold off
 plot(rrt256, rr256','+r'); hold on;
 plot(hrv256(:,1),hrv256(:,2),'+--b');
 title('RR intervals')
 ylabel('RR (s)')
 xlabel('(Removed RR intervals are plotted in red)')
end

fprintf('Now we are going to calculate an LF/HF-ratio\n for a 5min sliding window with 90 percent overlap\n Press any key to continue\n');
pause;

% calculate LF/HF ratio for a sliding window
out = lfhf_sliding_win(hrv256,300,0.9); %0.9 means take a 5min segment every 30 seconds!

% and plot it all
figure(2); hold off;
subplot(3,1,1)
plot(out.t/3600,out.lfhf)
axis([min(out.t/3600) max(out.t/3600) min(out.lfhf) max(out.lfhf)])
ylabel('lf/hf')
% title of graph
filename = 'LF/HF, StdDev and HR for 5 min window with 270s overlap';
title(filename)
subplot(3,1,2)
plot(out.t/3600,out.sdnn)
axis([min(out.t/3600) max(out.t/3600) min(out.sdnn) max(out.sdnn)])
ylabel('SD')
subplot(3,1,3)
plot(out.t/3600,60./out.hr)
axis([min(out.t/3600) max(out.t/3600) min(60./out.hr) max(60./out.hr)])
ylabel('HR')
xlabel('time (hours)');

fprintf('now we will look at higher sample rates and different thresholds\n');

%1024 Hz data
ecgInterp1024filt = interpft(ecg128filt,length(ecg128filt)*8);
% normal threshold
[rrs256 RS] = rrextract(ecgInterp1024filt',samplerate0*8,0.2,testmode);
fprintf('Notice the beat at 236.3 seconds has not been found,\n because the no. of samples in the integration window does not scale exactly\n Press any key to contniue\n');
pause;
[rrs256 RS] = rrextract(ecgInterp1024filt',samplerate0*8,0.1,testmode);
fprintf('Lowering the threshold from 0.2 to 0.1 allows us to catch this beat\n'); 


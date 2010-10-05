function plotWiener(model, noise, observation, utilda, Fs, Np)
% plotWiener(model, noise, observation, utilda, Fs)

if nargin<6
Np=4;
end

[Pm Fm] = psd(model,length(model),Fs,hamming(length(model)));
[Pn Fn] = psd(noise,length(noise),Fs,hamming(length(noise)));
[Psig Fsig] = psd(observation,length(observation),Fs,hamming(length(observation)));
[Pw Fw] = psd(utilda,length(utilda),Fs,hamming(length(utilda)));

figure; 
t=[1/Fs:1/Fs:length(observation)/Fs];
mn=min(min([model noise observation utilda]));  mn = mn - (0.1*mn);
mx=max(max([model noise observation utilda]));  mx = mx + (0.1*mx);

subplot(Np,1,1);plot(t,observation); ylabel('Unfiltered signal'); 
title('Wiener Filter'); 
subplot(Np,1,2);plot(t,model); axis([0 max(t) mn mx]); ylabel('Model');
subplot(Np,1,3);plot(t,noise); axis([0 max(t) mn mx]);ylabel('Noise');
subplot(Np,1,4);plot(t,utilda);axis([0 max(t) mn mx]);ylabel('Filtered Data');
if Np > 5
subplot(Np,1,5);plot(t,observation);hold on;plot(t,utilda,'r');axis([0 max(t) mn mx]);ylabel('Comparison1');
end
xlabel('time (s)')
figure;

%subplot(2,1,1);
%[Pm Fm] = psd(model,length(model),Fs,hamming(length(model)));
%[Pn Fn] = psd(noise,length(noise),Fs,hamming(length(noise)));
%[Psig Fsig] = psd(observation,length(observation),Fs,hamming(length(observation)));

%plot(Fm,10*log10(Pm),'g')
%hold on
%plot(Fn,10*log10(Pn),'m--')
%plot(Fs,10*log10(Ps),'k-.')
%xlabel('Frequency (Hz)')
%ylabel('Power (dB)')
%legend('Model','Noise','Unfiltered Signal') 
%axis([0 128 -80 19]);

%subplot(2,1,2);
%[Pm Fm] = psd(model,length(model),Fs,hamming(length(model)));
%[Pn Fn] = psd(noise,length(noise),Fs,hamming(length(noise)));
%[Ps Fs] = psd(observation,length(observation),Fs,hamming(length(observation)));
%[Pw Fw] = psd(utilda,length(utilda),Fs,hamming(length(utilda)));

%[Pw,Fw] = pwelch(utilda,hamming(length(wfiltest)),[], length(wfiltest),Fs); 

plot(Fm,10*log10(Pm),'g')
hold on
plot(Fn,10*log10(Pn),'m--')
plot(Fsig,10*log10(Psig),'k-.')
plot(Fw,10*log10(Pw),'r:')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend('Model','Noise','Unfiltered Signal','Filtered Signal') 
axis([0 128 -100 19]);



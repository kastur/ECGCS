function [yhat H] = wienerFilter(ideal,observation,R,graphicsFlagOn,Fs);
% 
% filtdata = wienerFilter(ideal,observation);
%
% FFT based Wiener filter in one dimension
% 
% Given a ideal of our perfect underlying signal that
% we wish to recover, we estimate the noise from
% noise = observation-ideal;
% The filtering is then performed in the frequency 
% domain by constructing the optimal (Wiener) filter
% for this noise/ideal estimate ... H=Sf2./(Sf2/Nf2)
% See Ref: Numerical Recipes in C, chapter 13.3, Press 
% Ref: http://www.phys.uni.torun.pl/nrbook/c13-3.pdf 
%
% G.D. Clifford 2004 gari AT mit DOT edu

% Further info:
% R=real(fft(r)) where r is the smearing function
% in the time domain. Default==1 ... no smearing.
% To prevent a graphical display of the results, 
% pass an extra argument==0
% For plotting, a sampling freq Fs=256 is assumed.
% You may pass a fifth arguement to change this
% (although it must be integer
%
% Licensing: GNU GPL applies


if nargin < 5
Fs=256;
end
Fs=round(Fs);

if nargin < 4
graphicsFlagOn=1;
end

if nargin < 3
R=1;     % response function     % ... no smearing
end

% estimate noise from ideal
noise = observation-ideal;


% work out how long to make FFT
N=length(observation);

% Wiener filter
Sf2=real(fft(ideal,N*2-1)).^2;   % Smeared ideal
Nf2=real(fft(noise,N*2-1)).^2;   % noise
Cf=real(fft(observation,N*2-1)); % ~= sqrt(Sf2+Nf2); % Corrupted ideal
      
H=Sf2./(Sf2+Nf2);              % Optimal filter

Yhat=H.*Cf/R;                  % 'uncorrupted' ideal estimate ...

yhat=real(ifft(Yhat));           % ..... in time domain

% ...compensate for FFT being two sided in matlab   
yhat=yhat(1:length(observation)); 

if graphicsFlagOn==1

plotWiener(ideal, noise, observation, yhat, Fs)

end

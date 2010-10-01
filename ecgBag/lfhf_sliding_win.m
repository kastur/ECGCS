function [out] = lfhf_sliding_win(rr_in, win, olap)

% [out] = lfhf_sliding_win(rr_in, win, olap)
% out is a structure such that:
% out.lfhf = lfhf ratio for win seconds, overlapped by factor 'olap' (>0 <1)
% and negative values indicate skipping segments between windows.
% Each window has a mean time point of out.t, with a number of points 
% equal to out.N
% out.hr is the mean RR interval in each window and out.sdnn is the standard 
% deviation of RR intervals in each window. out.respMag is the largest 
% frequency component between 0.15 and 0.4 Hz and out.resp is the 
% corresponding freuency.
%
% These routines are made available under the GNU general public license. 
% If you have not received a copy of this license, please download from 
% http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% G. Clifford : gari@ieee.org

if nargin < 3
   olap = 0.5;
end
if nargin < 2
   win = 300;
end 
if nargin < 1
   error('I need RR intervals!')
end

[a b] = size(rr_in);
if (a<2 | b<2)
     rri = rr_in;
     rrt = cumsum(rr_in);  
else
    rrt = rr_in(:,1);
    rri = rr_in(:,2);
end

% check to see if we have silly values 
olap = 1-olap;
if(olap>=1)
 error('olap should be _between_ 0 and 1');
end
% for Lomb periodogram
freq_vect =  [1/1024 : 1/1024 : 512 * 1/1024];

% find the start and end of each window
ind=1;
offset = rrt(1);
tstart(1)=1;
tend(1) = min( find(rrt> (offset+win+(win*olap*(ind-1))))  ) ;
while(offset+win+(win*olap*(ind))< rrt(length(rrt)) )
     tstart(ind+1)=min( find(rrt> offset+((win)*olap*(ind)))  ) ;
     tend(ind+1) = min( find(rrt> offset+(win+(win*olap*(ind))))  );
 ind=ind+1;
end

% loop over each window
for(ind=1:length(tend))

% intialise
 lfhf(ind)= NaN;
 t(ind)   = NaN;
 N(ind)   = NaN;
 hr(ind)  = NaN;
 sdnn(ind)= NaN;
 %resp(ind)= NaN;
 %respMag(ind)= NaN;

% define data window
 data = rri(tstart(ind):tend(ind)); % RR interval
 datat = rrt(tstart(ind):tend(ind));% time of RR intervals

% time at middle of each window
 t(ind)= mean([datat(1) datat(length(datat))]);
% average heart rate
 hr(ind)=mean(data);
% no. of points in each window
 N(ind)=length(data);

% have we got enough data/high enough HR to reach 0.4Hz
 if(mean(data)>(60/48))
     fprintf('warning - Peak Nyquist frequency not achieved up to 0.4Hz\n');
 end
 if((N/win)<0.8)
     fprintf('warning - Average Nyquist frequency not achieved up to 0.4Hz\n');
 end

 if(N<200 & N>=100)
     fprintf('warning - only %i points in window\n',N);
 end
 if(N<100)
     fprintf('warning - only %i points in window - statistically inaccurate results\n',N);
 end

 % SD (not strictly the SDNN!
 sdnn(ind) = std(data);

% PSD using unevenly sampled technique
 [Px, Prob] = lomb(datat, data-mean(data), freq_vect);
       
% check the statistical accuracy of the data       
 probFlag(ind)=1;
 if(isempty(Prob<0.05));
     fprintf('warning - no statistically relevant period fround in window\n');
     probFlag(ind)=0;
 end

% calculate the LF/HF ratio and the LF and HF (not normalised values)
 [lfhf(ind) lf(ind) hf(ind)] =  calc_lfhf(freq_vect,Px);

 fprintf('t=%f, lfhf=%f, sdnn=%f, hr=%f, N=%i \n',t(ind),lfhf(ind),sdnn(ind),hr(ind),N(ind));
end

% median filter the data with enough data to span beyond the overlap region
mfwin = round((1/(1-olap))+1);
medfiltlfhf = medfilt1(lfhf,mfwin);

% create an output structure of all our calculations
out.medfiltlfhf = medfiltlfhf;
out.lfhf = lfhf;
out.sdnn = sdnn;
out.t    = t;
out.hr   = hr;
out.N    = N;
out.lf = lf;
out.hf = hf;






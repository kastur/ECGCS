function [lfhf, lf, hf] =  calc_lfhf(Fx,Px)

% [lfhf] =  calc_lfhf(Fx,Px);
% Calculates a the LF/HF-ratio for a given (linear) PSD Px over 
% a given linear frequency range Fx
% Also:
% [lfhf lf hf] =  calc_lfhf(Fx,Px);
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


format long e

% check if defaults are changed
if nargin < 2
  error('frequency and power required: [lfhf lf hf] =  calc_lfhf(Fx,Px)')
end

%%%%%%%%%% set general defaults as per Malik-96 et al
LF_lo = 0.04;
LF_hi = 0.15;
HF_lo = 0.15; 
HF_hi = 0.4;
LF = 0.095;
HF = 0.275;
pow=0;
pow_cent=0;
peak_A=0;
peak_F=0;
j=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bin size of the PSD so that we can calculate the LF and HF metrics
binsize=Fx(2)-Fx(1);

% find the indexes corresponding to the LF and HF regions
indl = find( (Fx>=LF_lo) & (Fx<=LF_hi) );
indh = find( (Fx>=HF_lo) & (Fx<=HF_hi) );

% calculate metrics
lf   = binsize*sum(Px(indl));
hf   = binsize*sum(Px(indh));
lfhf =lf/hf;

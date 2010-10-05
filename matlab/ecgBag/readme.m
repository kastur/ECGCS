% ecgBag.tar is a collection of ECG processing algorithms for use
% with Matlab mostly by Gari Clifford and available under the GNU
% public license (see bottom of this file).  
% 
%
% They are not intended as optimal algorithms, but they're not bad.....
% Please don't blame me if anything goes wrong - but feel free to
% email me and ask questions, which I'll answer if I have time
% Gari Clifford (gari@ieee.org) 
%
% calc_lfhf.m   - calculates the LF/HF ratio, LF and HF metrics       
% lomb.m        - calculates the Lomb-Scargle periodogram 
%                 (PSD for unevenly sampled signals)  
% rpeakdetect.m - finds R peaks in ECG - see rrextract
% rrextract.m   - runs over rpeakdetect to extract RR intervals
% clean_hrv4.m  - removes abnormally timed beats     
% clean_RR_RS.m - removes abnormally timed & shaped beats     
% example.m     - an example file to run all this     
% maxjiggle.m   - finds peaks if the fiducial point is slightly off
% interp_RR.m   - interpolates an unevenly sampled signal by your favourite 
%                 (ad hoc) method. NOT recommended - it just introduces 
%                 errors into a PSD calculation - use the Lomb-Scargle 
%                 method, it's so much better.
% filterECG256Hz.m    - ECG, 256Hz, Least Sq. FIR LP & HP filters 
%                       (cascaded) 70dB 0.05-40Hz 1dB ripple
% filterECG128Hz.m    -  ECG, 128Hz, Least Sq. FIR LP & HP filters 
%                       (cascaded)  70dB 0.05-40Hz 1dB ripple
% lfhf_sliding_win.m  - calculates a time series of LF/HF values by 
%                       sliding a window across the data
%                      sliding a window across the data
%
% generate_rr_with_FM.m       - Creates an artificial RR interval
%                               time series
%
% down_sample_realistically.m - takes an evenly sampled time series
%                               and extracts a plausible (irregularly
%                               sampled) RR interval time series
% ECGfilterTimeAlignedSVD.m   - uses a QRS detector and SVD to filter ECG
% select_train.m              - makes a matrix of fiducual point centered
%                               P-QRS-T complexes.
% recon_ECG.m                 - unwraps a matrix of fiducual point centered
%                               P-QRS-T complexes into a 1-D signal
% svdFilter.m                 - filters data using SVD - variant of ECGfilterTimeAlignedSVD.m
% parabolic_filter.m          - Parabolic Filter for baseline subtraction
% wienerFilter.m              - Wiener filtering for the ECG
% plotWiener.m                - plots the spectral response of wienerFilter.m
% wavelet_decompose.m         - A wavelet filter for ECG analysis by Nick Hughes
% wabp.m                      - James Sun's version of Wei Zong's wabp.c (on www.physionet.org)
%                                   (i.e. - a QRS detector for blood pressure signals)
%
% snr.m                       - calcualtes the signal to noise ration in dB
% 
% conv_norm_annot.c  ... this C-prog parses the output of the WFDB 
% rdann program to extract the Normal to Normal intervals and puts 
% it in a readable format for Matlab - read the header for details
% conv_sleep_annot.c  ... as above, but works on sleep staging...  
% - again, read the header for details, and type:
% type readme 
%
% This software is offered freely and without warranty under the 
% GNU public license. Please don't blame me if anything goes wrong :)
%
% gari@ieee.org
 
% for the polysomnographic DB, you would run rdann as
% rdann -r slp02a -a st
% and it spits out the data to stdout .. . so best to do -
% rdann -r slp02a -a st > fileout.dat
% (note that you must have the data on a current path to do this).
%
% put simply:
% rdann -r slp01a -a st > slp01a.ascii
% conv_sleep_annot  slp01a.ascii  slp01a.tRnKart  > log1
% rdann -r slp01a -a ecg > slp01a.tmp
% conv_norm_annot slp01a.tmp slp01a.bts > log2
%
% slp01a.tRnKart should be a 3 column file of [time_stamp sleep_stage other_label]
% slp01a.bts  should be a 2 column file of [time_stamp beat_label]
%
% The time stamps are in seconds since the beginning of the file 
% (so check your sampling frequency is the same as the #define in the C file).
%
% The sleep stage are:
%    -1  W   subject is awake
%     0  R   REM sleep
%     1  1   sleep stage 1
%     2  2   sleep stage 2
%     2  2   sleep stage 3
%     4  4   sleep stage 4
%
% The other labels are
%  -0.1     No label in 3rd column
%    -2  H   Hypopnea
%  -0.2  HA  Hypopnea with arousal
%    -3  OA  Obstructive apnea
%  -0.3  X   Obstructive apnea with arousal
%    -4  CA  Central apnea
%  -0.4  CAA Central apnea with arousal
%    -5  L   Leg movements
%  -0.5  LA  Leg movements with arousal
%  -0.6  A   Unspecified arousal
%   -10  MT  Movement time
%
% The beat labels are
% 1 == normal
% 2 -> 19 == abnormal/ectopics 
% 2 = L
% 3 = R 
% 4 = B 
% 5 = A 
% 6 = a 
% 7 = J 
% 8 = S 
% 9 = V 
% 10 = r 
% 11 = F 
% 12 = e 
% 13 = j 
% 14 = n 
% 15 = E 
% 16 = f 
% 17 = Q 
% 30 and above == artefact
% 20 -- unspecified (\?)
% 25 -- paced (P)
% 0 -- error - classification not recognized

%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 
% Most of these m-files and their dependents should be freely available from 
% Physionet -  http://www.physionet.org/ - in the near future... but please
% report any bugs to gari@ieee.org



% Change Log
% ----------
% Robert Tratnig at Graz Uni added the following comments and changes:
%
% clean_hrv_4.m:
% I have added a line that says "stp=0;".
% I call this line "38".
% According to this indication I added lines 38,39,40,43,44,45,47,50,51,82.
% Additionally I have altered line 40 (I added a condition: "(stp==0) & ".
% 
% The same scheme I have applied to clean_RR_RS.m
% Here I call the line that says "stp=0;" line "51".
% In that context I added lines 51,52,53,56,57,58,60,61,62,102.
% Also here I have altered a line: 54 (I added a condition: "(stp==0) & ".
% 
% I have to do this because of the following reason:
% When analyzing the file cu06.dat of the Creighton University Databank (CU_DB)
% in the MIT-data collection, the parameter transfered to this function
% clean_hrv4.m and clean_RR_RS.m (temp_hrv) was just a 3 rows long array of 2
% columns. Additionally no row satisfied the condition of the while loop, so I
% had an error called "Index exceeds matrix dimensions". 
%
% Many thanks for that Rob!
%
%
%

% ToDo
% ----
%
% - vectorise lomb.m ... currently very slow. this code was emailed to
%   me sometime in the past by a colleague who obtained it freely from the
%   internet. I have made minor modifications but have so far failed to have 
%   time (or sufficient reason) to vectorise it (lomb.c does the job!)
%
% - Add optimisation routine to make either end of RR tachogram window 
%   have same mean and gradient within a tolerance as per Schreiber's
%   surrogate window matching. 









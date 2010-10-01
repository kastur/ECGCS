% WAVELET_DECOMPOSE  Decompose a signal upto a certain scale with a given wavelet
%
% [approximations, details] = wavelet_decompose(signal, scale, wavelet)
%
% signal - 1-D real valued input signal for decomposition
% scale - scale to filter upto
% wavelet - wavelet basis to use for decomposition (character string)
%
% approximations - matrix of successive approximation signals (column-wise)
% details - matrix of successive detail signals (column-wise)
%
% Available wavelet names 'wname' are:
%    Daubechies: 'db1' or 'haar', 'db2', ... ,'db45'
%    Coiflets  : 'coif1', ... ,  'coif5'
%    Symlets   : 'sym2' , ... ,  'sym8', ... ,'sym45'
%    Biorthogonal:
%        'bior1.1', 'bior1.3' , 'bior1.5'
%        'bior2.2', 'bior2.4' , 'bior2.6', 'bior2.8'
%        'bior3.1', 'bior3.3' , 'bior3.5', 'bior3.7'
%        'bior3.9', 'bior4.4' , 'bior5.5', 'bior6.8'.
%    Reverse Biorthogonal: 
%        'rbio1.1', 'rbio1.3' , 'rbio1.5'
%        'rbio2.2', 'rbio2.4' , 'rbio2.6', 'rbio2.8'
%        'rbio3.1', 'rbio3.3' , 'rbio3.5', 'rbio3.7'
%        'rbio3.9', 'rbio4.4' , 'rbio5.5', 'rbio6.8'.
%
% e.g. ...
% [approximations, details] = wavelet_decompose(signal, 5, 'bior3.3')
%
% The function returns an
% approximation matrix and a details matrix. The columns of the approx
% matrix are the filtered signals, so column 'i' has the filtered signal at
% level 'i' (which is "scale" 2^i). As the level increases so too does the
% degree of filtering. The columns of the details signal give the part of
% the signal that was removed in the filtering process at that level.
% Try different wavelets: bior4.4, bior5.5, bior6.8 are all good for ECG.
%
% Copyright (c) Nick Hughes (2001)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [approximations, details] = wavelet_decompose(signal, scale, wavelet)

sig_length = length(signal);

approximations = zeros(sig_length, scale);
details = zeros(sig_length, scale);

[C,L] = wavedec(signal, scale, wavelet);

for i=1:scale,
  approximations(:,i) = wrcoef('a', C, L, wavelet, i);
  details(:,i) = wrcoef('d', C, L, wavelet, i);
end


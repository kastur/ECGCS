function [m, i] = localmax2(x, w)
% LOCALMAX  Find indices and amplitudes of local maxima 
% [m,i] = localmax2(x, w) returns the indices and maxima defined 
% over a local window of size 2w+1 given by w points on either 
% side of the point being considered as a local maximum..
%
% P.E.McSharry
% These routines are made available under the GNU general public license. 
% If you have not received a copy of this license, please download from 
% http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% G. Clifford : gari@ieee.org

N = length(x);

k = 2*w+1;
y = zeros(k,1);

l = 0;
for j=w+1:N-w
   y = x(j-w:j+w);
   [ymax,imax] = max(y);
   if imax == w+1
      l = l+1;
      m(l) = ymax;
      i(l) = j; 
   end
end

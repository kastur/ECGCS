function imax = maxjiggle(x, ind, d)
% imax = maxjiggle(x, ind, d) returns the indices of the local maxima by 
% jiggling the initial indices ind by +- d. 
% by P.E.McSharry
% These routines are made available under the GNU general public license. 
% If you have not received a copy of this license, please download from 
% http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% G. Clifford : gari@ieee.org

N = length(x)
Ni = length(ind)

xmin = min(x);
X = ones(Ni,2*d+1)*xmin;

for i=-d:d
   j = find(1 <= ind+i & ind+i <= N);
   X(j,d+i+1) = x(ind(j)+i); 
end

[vmax,ivmax] = max(X,[],2);
imax = ind + ivmax-d-1;

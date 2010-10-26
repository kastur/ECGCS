/*
File Name: MIDWT.c
Last Modification Date:	06/14/95	13:01:15
Current Version: MIDWT.c	2.4
File Creation Date: Wed Oct 12 08:44:43 1994
Author: Markus Lang  <lang@jazz.rice.edu>

Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
Created by Markus Lang, Department of ECE, Rice University. 

This software is distributed and licensed to you on a non-exclusive 
basis, free-of-charge. Redistribution and use in source and binary forms, 
with or without modification, are permitted provided that the following 
conditions are met:

1. Redistribution of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer.
2. Redistribution in binary form must reproduce the above copyright notice, 
   this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software 
   must display the following acknowledgment: This product includes 
   software developed by Rice University, Houston, Texas and its contributors.
4. Neither the name of the University nor the names of its contributors 
   may be used to endorse or promote products derived from this software 
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY WILLIAM MARSH RICE UNIVERSITY, HOUSTON, TEXAS, 
AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE UNIVERSITY 
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTIONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE), PRODUCT LIABILITY, OR OTHERWISE ARISING IN ANY WAY OUT OF THE 
USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For information on commercial licenses, contact Rice University's Office of 
Technology Transfer at techtran@rice.edu or (713) 348-6173

Change History: Fixed the code such that 1D vectors passed to it can be in
                either passed as a row or column vector. Also took care of 
		the code such that it will compile with both under standard
		C compilers as well as for ANSI C compilers
		Jan Erik Odegard <odegard@ece.rice.edu> Wed Jun 14 1995

		Fix minor bug to allow maximum number of levels

decription of the matlab call:
%y = midwt(x,h,L);
% 
% function computes the inverse discrete wavelet transform y for a 1D or 2D
% input signal x.
%
%    Input:
%	x    : finite length 1D or 2D input signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%
% see also: mdwt, mrdwt, mirdwt

*/

#include <math.h>

#define max(A,B) (A > B ? A : B)
#define mat(a, i, j) a[m*j+i]  /* macro for matrix indices */

template<class T>
void bpsconv(vector<T>& x_out, int lx, vector<T>& g0, vector<T>& g1, int lhm1, 
	int lhhm1, vector<T>& x_inl, vector<T>& x_inh);

template<class T>
void MIDWT(const vector<T>& y, const vector<T>& h, int L, vector<T>& x) {
  int lh = h.size();
  int lhm1 = lh - 1;
  int lhhm1 = lh/2 - 1;
  int n = y.size();

  vector<T> xdummy(n);
	vector<T> ydummyl(n+lh/2-1);
	vector<T> ydummyh(n+lh/2-1);
	vector<T> g0(lh);
	vector<T> g1(lh);

	int c_o_a, ic;

	int actual_n;

  /* synthesis lowpass and highpass */
  for (int i = 0; i < lh; i++) {
    g0[i] = h[i];
    g1[i] = h[lh-i-1];
  }
  for (int i = 1; i <= lh; i += 2)
    g1[i] = -g1[i];
  

  /* 2^L */
  int sample_f = 1;
  for (int i = 1; i < L; i++)
    sample_f = sample_f*2;
  
  actual_n = n/sample_f;

  for (int i = 0; i < n; i++)
    x[i] = y[i];
  
  /* main loop */
  for (int actual_L = L; actual_L >= 1; actual_L--) {
    c_o_a = actual_n/2;
    /* store in dummy variable */
    ic = c_o_a;
    for  (int i = 0; i < c_o_a; i++) {
			ydummyl[i+lhhm1] = x[i];  
			ydummyh[i+lhhm1] = x[ic];
      ic++;
    }
    /* perform filtering lowpass and highpass*/
    bpsconv(xdummy, c_o_a, g0, g1, lhm1, lhhm1, ydummyl, ydummyh); 
    /* restore dummy variables in matrices */
    for (int i=0; i<actual_n; i++)
      x[i] = xdummy[i];  
    actual_n = actual_n*2;
  }
}

template <class T>
void bpsconv(vector<T>& x_out, int lx, vector<T>& g0, vector<T>& g1, int lhm1, 
	int lhhm1, vector<T>& x_inl, vector<T>& x_inh) {
  T x0, x1;
  for (int i = lhhm1 - 1; i > -1; i--) {
    x_inl[i] = x_inl[lx+i];
    x_inh[i] = x_inh[lx+i];
  }

  int tj;
  for (int i = 0, ind = 0; i < lx; i++) {
    x0 = 0;
    x1 = 0;
    tj = -2;
    for (int j=0; j <= lhhm1; j++){
      tj+=2;
      x0 = x0 + x_inl[i+j]*g0[lhm1-1-tj] + x_inh[i+j]*g1[lhm1-1-tj] ;
      x1 = x1 + x_inl[i+j]*g0[lhm1-tj] + x_inh[i+j]*g1[lhm1-tj] ;
    }
    x_out[ind++] = x0;
    x_out[ind++] = x1;
  }
}

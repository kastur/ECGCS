#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "compressed-sensing-inl.h"

float x[] = { 
   118,
   117,
   117,
   117,
   117,
   117,
   117,
   117,
   118,
   120,
   122,
   123 };

int main() {
	typedef float T;
	unsigned int lfsr_state = 854875398;
	
	const int N = 12;
  const int K = 10;
  float y[K];
  float xhat[N];
	encode_bernoulli(lfsr_state, N, K, x, y);

 
	{
		vector<T> temp(K, y);
		PrintT(temp);
	}

	decode_dwt(lfsr_state, K, N, y, xhat);
	{
		vector<T> temp(N, xhat);
		PrintT(temp);
	}


	return 0;
}

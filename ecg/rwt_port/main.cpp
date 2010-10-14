// Author: Kasturi Raghavan (kasturir@ucla.edu)
// \file
// Tests out dwt and idwt functions provided by the Rice
// DWT package. The wavelet scaling coefficients must be provided
#include <iostream>
#include <cstdlib>
#include "mdwt_r.h"
#include "midwt_r.h"

using namespace std;

int main() {
	const int m = 1;
	const int n = 4;
	const int scaling_coef_n = 4;
	const int levels_n = 1;
	const double scaling_coef[] = {0.4830, 0.8365, 0.2241, -0.1294};
	const double signal_t[] = { 1.0, 1.0, 1.0, 1.0 };
	
	double *signal_dwt = (double*)calloc(m*n, sizeof(double));
	double *signal_t1 = (double*)calloc(m*n, sizeof(double));

	MDWT<m, n, scaling_coef_n, levels_n>(signal_t, scaling_coef, signal_dwt);
	MIDWT<n, m, scaling_coef_n, levels_n>(signal_t1, scaling_coef, signal_dwt);

	cout << "Original signal: " << endl;
	for (int i = 0; i < m*n; i++)
		cout << signal_t[i] << endl;

	cout << endl;

	cout << "DWT of original signal: " << endl;
	for (int i = 0; i < m*n; i++)
		cout << signal_dwt[i] << endl;

	cout << endl;

	cout << "IDWT, recovered original signal: " << endl;
	for (int i = 0; i < m*n; i++)
		cout << signal_t1[i] << endl;
}

#include <iostream>
#include <cstdlib>
#include "mdwt_r.h"
#include "midwt_r.h"

using namespace std;

int main() {


	const int m = 4;
	const int n = 1;
	double scaling_coef[] = {0.4830, 0.8365, 0.2241, -0.1294};
	const int scaling_coef_n = 4;
	const int levels_n = 1;

	double *signal_dwt = (double*)calloc(m*n, sizeof(double));
	double *signal_t1 = (double*)calloc(m*n, sizeof(double));
	double signal_t[] = { 1.0, 1.0, 1.0, 1.0 };

	MDWT(signal_t, m, n, scaling_coef, scaling_coef_n, levels_n, signal_dwt);

	MIDWT(signal_t1, m, n, scaling_coef, scaling_coef_n, levels_n, signal_dwt);


	for(int i = 0; i < m*n; i++)
		cout << signal_dwt[i] << endl;

	cout << endl;

	for(int i = 0; i < m*n; i++)
		cout << signal_t1[i] << endl;

}

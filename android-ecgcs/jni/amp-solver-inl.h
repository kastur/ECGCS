#ifndef AMP_SOLVER_INL_H
#define AMP_SOLVER_INL_H

#include <math.h>
#include "linalg-inl.h"
#include "largest-element-inl.h"

using namespace linalg;

namespace amp_solver {


template <class T>
T eta(T x, T threshold) {
  if (x > threshold)
    return x - threshold;
  else if (x < -threshold)
    return x + threshold;
  else
   	return 0;
}

template <class T>
T eta_prime(T x, T threshold) {
  if ((x > threshold) || (x < -threshold))
    return 1;
  else
    return 0;
}

template <class T>
void Eta(const vector<T>& x, T threshold, vector<T>& y) {
  for (size_t i = 0; i < y.size(); i++)
    y[i] = eta(x[i], threshold);
}

template <class T>
void EtaPrime(const vector<T>& x, T threshold, vector<T>& y) {
  for (size_t i = 0; i < y.size(); i++)
    y[i] = eta_prime(x[i], threshold);
}

template<class T, class M>
void solve(const M& A, const vector<T>& y, size_t n_iterations, T tol,
  vector<T>& xhat) {
  const T one = 1;
  size_t K = y.size();
  size_t N = A.size() / K;
 // y = A x  A is k x n
  vector<T> z(K);
	T y_norm;

	// DEBUG {
	vector<T> norms(n_iterations);
  // }

	vector<T> Atz(N);  // A' * z
  vector<T> gamma(N);
  vector<T> abs_gamma(N);
  T threshold;
	vector<T> Axhat(K);  // A * xhat
  vector<T> resid(K);  // A * xhat - y
	T resid_norm;
  vector<T> threshold_gamma(N);
	T sum_threshold_gamma;
  vector<T> weighted_z(K);
  
  // Initial values
  // Zeros(xhat);  // should we use some other initial guess?
	// WARNING: caller initialized xhat initial value!
  Copy(y, z);
	Norm2(y, y_norm);

	for (size_t t = 0; t < n_iterations; t++) {
    TimesTranspose(A, z, Atz);
		Plus(xhat, Atz, gamma);
		Abs(gamma, abs_gamma);
		largestElement(abs_gamma, K, threshold);
		Eta(gamma, threshold, xhat);
		Times(A, xhat, one, Axhat);
		Plus(y, Axhat, one, -one, resid);

		Norm2(resid, resid_norm);

		EtaPrime(gamma, threshold, threshold_gamma);
		Sum(threshold_gamma, sum_threshold_gamma);
		Times(z, (sum_threshold_gamma/K), weighted_z);
		Plus(resid, weighted_z, z);

		norms[t] = resid_norm/y_norm;
	}

	// PrintT(norms);
}


}  // end namespace amp_solver

#endif  // AMP_SOLVER_INL_H

#ifndef COMPRESSED_SENSING_INL_H_
#define COMPRESSED_SENSING_INL_H_

#include "linalg-inl.h"
#include "linalg-debug-inl.h"
#include "lfsr-inl.h"
#include "function-matrix-inl.h"
#include "amp-solver-inl.h"

using amp_solver::solve;

template<class T>
int encode_bernoulli(unsigned int lfsr_state, size_t N, size_t K, T* x_data, T* y_data) {
	LFSRRand rand(lfsr_state);
	DwtFunctionMatrix<T> A(rand, N, K);
	vector<T> x(N, x_data);
	vector<T> y(K);
	Times(A.GetMatrix(), x, y);
	memcpy(y_data, &(y[0]), (sizeof(T)*K));
	return rand.getState();
}

template<class T>
int decode_dwt(unsigned int lfsr_state, size_t K, size_t N, T* y_data, T* xhat_data, size_t iter = 100, T tol = 0.01) {
	LFSRRand rand(lfsr_state);
  DwtFunctionMatrix<T> A(rand, N, K);

  vector<T> y(K, y_data);
  vector<T> zhat(N);
	zhat.fill(0);

  solve(A, y, iter, tol, zhat);

	vector<T> xhat(N);
	A.ApplyFunction(zhat, xhat);

	memcpy(xhat_data, &(xhat[0]), (sizeof(T)*N));
  return rand.getState();
}

#endif  // COMPRESSED_SENSING_INL_H_

#ifndef DECODE_INL_H_
#define DECODE_INL_H_

#include "linalg-inl.h"
#include "linalg-debug-inl.h"
#include "largest-element-inl.h"
#include "lfsr-inl.h"
#include "amp-solver-inl.h"
#include "midwt_r.h"
#include "mdwt_r.h"


using amp_solver::solve;

const float dwt_coef_arr[4] = {
0.48296291314453415610685738101893,
0.83651630373780794247551284570363,
0.22414386804201338887487793272157,
-0.12940952255126039749377753196313 };

const vector<float> dwt_coef(4, dwt_coef_arr);
const int dwt_levels = 1;

template <class T>
class FunctionMatrix {
 public:
  FunctionMatrix(LFSRRand& rand, int m, int n) : matrix_(m*n) {
		T norm_factor = m;
    BernoulliNormalized(rand, norm_factor, matrix_);
  }

  void ApplyFunction(const vector<T>& x, vector<T>& y) const {
    y.fill(0);
    MIDWT<T>(x, dwt_coef, dwt_levels, y);
  }

  void ApplyTransposeFunction(const vector<T>& x, vector<T>& y) const {
    y.fill(0);
    MDWT<T>(x, dwt_coef, dwt_levels, y);
  }

  const vector<T>& GetMatrix() const {
    return matrix_;
  }
  
  size_t size() const {
    return matrix_.size();
  }

 private:
  vector<T> matrix_;
};

template <class T>
void Times(const FunctionMatrix<T>& R, const vector<T>& x, vector<T>& z) {
  vector<T> y(x.size());
	R.ApplyFunction(x, y);
  Times(R.GetMatrix(), y, z);
}

template <class T>
void TimesTranspose(const FunctionMatrix<T>& R, const vector<T>& x, vector<T>& z) {
  vector<T> y(R.GetMatrix().size() / x.size());
  TimesTranspose(R.GetMatrix(), x, y);
	R.ApplyTransposeFunction(y, z);
}

template <class T>
void Times(const FunctionMatrix<T>& R, const vector<T>& x, T a, vector<T>& z) {
  vector<T> y(x.size());
	R.ApplyFunction(x, y);
  Times(R.GetMatrix(), y, a, z);
}

template <class T>
void Print(const FunctionMatrix<T>& F, int K) {
  Print(F.GetMatrix(), K);
}

template <class T>
void PrintT(const FunctionMatrix<T>& F, int K) {
  PrintT(F.GetMatrix(), K);
}

template<class T>
int encode(unsigned int lfsr_state, size_t N, size_t K, T* x_data, T* y_data) {
	LFSRRand rand(lfsr_state);
	FunctionMatrix<T> A(rand, N, K);
	vector<T> x(N, x_data);
	vector<T> y(K);
	Times(A.GetMatrix(), x, y);
	memcpy(y_data, &(y[0]), (sizeof(T)*K));
	return rand.getState();
}

template<class T>
int decode(unsigned int lfsr_state, size_t K, size_t N, T* y_data, T* xhat_data, size_t iter = 100, T tol = 0.01) {
	LFSRRand rand(lfsr_state);
  FunctionMatrix<T> A(rand, N, K);

  vector<T> y(K, y_data);
  vector<T> zhat(N);
	zhat.fill(0);

  solve(A, y, iter, tol, zhat);

	vector<T> xhat(N);
	A.ApplyFunction(zhat, xhat);

	memcpy(xhat_data, &(xhat[0]), (sizeof(T)*N));
  return rand.getState();
}

#endif  // DECODE_INL_H_

#ifndef FUNCTION_MATRIX_INL_H_
#define FUNCTION_MATRIX_INL_H_

#include "linalg-inl.h"
#include "linalg-debug-inl.h"
#include "lfsr-inl.h"
#include "midwt_r.h"
#include "mdwt_r.h"

// Daub. coefficients, see Rice Wavelet Toolbox
const float dwt_coef_arr[4] = {
0.48296291314453415610685738101893,
0.83651630373780794247551284570363,
0.22414386804201338887487793272157,
-0.12940952255126039749377753196313 };
const vector<float> dwt_coef(4, dwt_coef_arr);
const int dwt_levels = 1;

// DwtFunctionMatrix serves is an efficient implementation of
//   A = [SketchMatrix * InverseDwt * Identity] * x
//   A' = [SketchMatrix * InverseDwt * Identity]' * x
// Actually implemented with DWT and IDWT as a functions:
//   A = SketchMatrix * InverseDwt(x)
//   A' = Dwt(x) * SketchMatrix'
template <class T>
class DwtFunctionMatrix {
 public:
  DwtFunctionMatrix(LFSRRand& rand, int m, int n) : matrix_(m*n) {
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

namespace linalg {

template <class T>
void Times(const DwtFunctionMatrix<T>& R, const vector<T>& x, vector<T>& z) {
  vector<T> y(x.size());
	R.ApplyFunction(x, y);
  Times(R.GetMatrix(), y, z);
}

template <class T>
void TimesTranspose(const DwtFunctionMatrix<T>& R, const vector<T>& x, vector<T>& z) {
  vector<T> y(R.GetMatrix().size() / x.size());
  TimesTranspose(R.GetMatrix(), x, y);
	R.ApplyTransposeFunction(y, z);
}

template <class T>
void Times(const DwtFunctionMatrix<T>& R, const vector<T>& x, T a, vector<T>& z) {
  vector<T> y(x.size());
	R.ApplyFunction(x, y);
  Times(R.GetMatrix(), y, a, z);
}

}

template <class T>
void Print(const DwtFunctionMatrix<T>& F, int K) {
  Print(F.GetMatrix(), K);
}

template <class T>
void PrintT(const DwtFunctionMatrix<T>& F, int K) {
  PrintT(F.GetMatrix(), K);
}

#endif  // FUNCTION_MATRIX_INL_H_

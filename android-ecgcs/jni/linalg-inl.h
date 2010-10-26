#ifndef LINALG_INL_H_
#define LINALG_INL_H_

#include <math.h>
#include "vector-inl.h"

namespace linalg {

template <class T>
void Copy(const vector<T>& x, vector<T>& y) {
	y = x;
}

template <class T>
void Zeros(vector<T>& x) {
	x.fill(0);
}

template <class T>
void Eye(size_t n, vector<T>& A) {
	Zeros(A);
  for(size_t i = 0; i < n; i++)
    A[i*n + i] = 1;
}

template <class T>
void Times(const vector<T>& x, T a, vector<T>& y) {
  for (size_t i = 0; i < x.size(); i++) {
    y[i] = x[i] * a;
  }
}

template <class T>
void Times(const vector<T>& A, const vector<T>& x, vector<T>& y) {
  size_t n = x.size();
  size_t m = A.size() / n;
  for (size_t i = 0; i < m; i++) {
    T dot_product = 0;
    for (size_t j = 0; j < n; j++) {
      dot_product += A[i*n + j] * x[j];
    }
    y[i] = dot_product;
  }
}

template <class T>
void TimesTranspose(const vector<T>& A, const vector<T>& x, vector<T>& y) {
  size_t n = x.size();
  size_t m = A.size() / n;
  for (size_t j = 0; j < m; j++) {
    T dot_product = 0;
    for (size_t i = 0; i < n; i++) {
      dot_product += A[i*m + j] * x[i];
    }
    y[j] = dot_product;
  }
}

template <class T>
void Times(const vector<T>& A, const vector<T>& x, T a, vector<T>& y) {
  size_t n = x.size();
  size_t m = A.size() / n;
  for (size_t i = 0; i < m; i++) {
    T dot_product = 0;
    for (size_t j = 0; j < n; j++) {
      dot_product += A[i*n + j] * a * x[j];
    }
    y[i] = dot_product;
  }
}

template <class T>
void Plus(const vector<T>& x, const vector<T>& y, T a, T b, vector<T>& z) {
	for (size_t i = 0; i < x.size(); i++) {
		z[i] = a*x[i] + b*y[i];
	}
}

template <class T>
void Plus(const vector<T>& x, const vector<T>& y, vector<T>& z) {
	for (size_t i = 0; i < x.size(); i++) {
		z[i] = x[i] + y[i];
	}
}

template <class T>
void Plus(const vector<T>& x, T y, vector<T>& z) {
	for (size_t i = 0; i < x.size(); i++) {
		z[i] = x[i] + y;
	}
}

template <class T>
void Sum(const vector<T>& x, T& sum) {
	sum = 0;
	for (size_t i = 0; i < x.size(); i++)
		sum += x[i];
}

template <class T>
void Abs(const vector<T>& x, vector<T>& y) {
  for (size_t i = 0; i < x.size(); i++)
    y[i] = fabs(x[i]);
}

template <class T>
void Norm2(const vector<T>& x, T& norm2) {
	T ssum = 0;
	for (size_t i = 0; i < x.size(); i++) {
		ssum += x[i]*x[i];
	}
	norm2 = sqrt(ssum);
}

}  // end namespace linalg 

#endif  // LINALG_INL_H_

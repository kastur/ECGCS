#ifndef LINALG_DEBUG_INL_H
#define LINALG_DEBUG_INL_H
#include <stdlib.h>
#include <stdio.h>

#define DEBUG 1

#ifdef DEBUG
#  ifdef ANDROID
#    include <android/log.h>
#    define D(x...) __android_log_print(ANDROID_LOG_DEBUG, "ECGCS", x)
#  else
#    define D(x...) printf(x)
#  endif
#else
#  define D(...) do { } while(0)
#endif

namespace linalg {

void Print(const char* str) {
  D("%s\n", str);
}

// Print to stdout
void Print(float x) {
	D("%+0.3f  ", x);
}

// Print to stdout m x n, m rows
template <class T>
void Print(const vector<T>& A, size_t m = 1) {
  size_t n = A.size() / m;
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
  		D("%+0.3f  ", A[i*n + j]);
    }
		D("\n");
  }
	D("\n");
}

// Print to stdout m x n transposed, i.e. n rows
template <class T>
void PrintT(const vector<T>& A, size_t m = 1) {
  size_t n = A.size() / m;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
  		D("%+0.3f  ", A[j*n + i]);
    }
		D("\n");
  }
	D("\n");
}


}  // end namespace linalg

#endif  // LINALG_DEBUG_INL_H

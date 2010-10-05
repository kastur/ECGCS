#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace std;

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
T nth_largest(const vector<T>& signal, size_t k) {
  vector<T> max_heap = signal;
  make_heap(max_heap.begin(), max_heap.end());
  size_t n = signal.size();
  T val = max_heap.front();
  if (n > k) {
    size_t j = n - k;
    for(size_t x = 0; x < j; x++) {
      pop_heap(max_heap.begin(), max_heap.end());
      max_heap.pop_back();
      val = max_heap.front();
    }
  }
  return val;
}

int test() {
  vector<int> signal;
  for(int i = -10; i <= 10; i++) {
    int v = eta(i, 4);
    signal.push_back(v);
    cout << eta(i, 4) << ", " << eta_prime(i, 4) << endl;
  }
  for(int i = signal.size(); i > 0; i--)
    cout << i << " largest " << nth_largest(signal, i) << endl;
  return 0;
}

template <class T>
void Ax(const vector<T>& A, const vector<T>& x, vector<T>* y) {
  int n = x.size();
  int m = A.size() / n;
  y->clear();
  y->resize(m, 0);
  for (int i = 0; i < m; i++) {
    T dot_product = 0;
    for (int j = 0; j < n; j++) {
      if (A[i*n + j] == 0) continue;
      dot_product += A[i*n + j] * x[j];
    }
    y->at(i) = dot_product;
  }
}

template <class T>
void Plus(const vector<T>& a, const vector<T>& b, vector<T>* c) {
  c->clear();
  c->resize(a.size(), 0);
  for (int i = 0; i < a.size(); i++) {
    c->at(i) = a[i] + b[i];
  }
}

template <class T>
void eye(int n, vector<T>* A) {
  A->clear();
  A->resize((n*n), 0);
  for(int i = 0; i < n; i++)
    A->at(i*n + i) = 1;
}

template <class T>
void zeros(int n, vector<T>* x) {
 x->clear();
 x->resize(n, 0);
}

template <class T>
void Abs(vector<T>* x) {
  for (int i = 0; i < x->size(); i++)
    x->at(i) = abs(x->at(i));
}

template <class T>
void Eta(vector<T>* x, T threshold) {
  for (int i = 0; i < x->size(); i++)
    x->at(i) = eta(x->at(i), threshold);
}

template <class T>
void rand_vec(int n, vector<T>* x) {
  x->clear();
  x->resize(n, 0);
  for (int i = 0; i < n; i++)
    x->at(i) = rand();
}

template <class T>
void print_vec(const vector<T>& x) {
  for (int i = 0; i < x.size(); i++)
    cout << x[i] << " ";
  cout << endl;
}

template <class T>
void print_mat(const vector<T>& A, size_t n) {
  int m = A.size() / n;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      cout << A[i*n + j] << ",";
    }
    cout << endl;
  }
}

template <class T>
void read_mat(ifstream& value_stream, int m, int n, vector<T>* A) {
  A->resize((m*n), 0);
  float val;
  char c;
  for(int i = 0; i < (m*n); i++) {
   value_stream >> val >> c;
   A->at(i) = val;
  }
}

template<class T>
void amp(const vector<T>& A, const vector<T>& y) {

  int K = y.size();
  int N = A.size() / K;

  int Q = 500; // mx iterations
  float tol = 0.0001; // tolerance

  vector<T> xhat;  // solution
  vector<T> z;
  vector<T> gamma;
  vector<T> temp;



// initial estimate
  zeros(N, &xhat);
  z = y;

//for (int t = 0; t < T; t++) {
  Ax(A, z, &temp);
  print_mat(A, N); 
  Plus(xhat, temp, &gamma);
  Abs(&gamma);
  T threshold = nth_largest(gamma, K);
  Eta(&gamma, threshold);
//}
}

template <class T>
void init() {
}

int main() {
  typedef float T;
  int K = 40;
  int N = 64;

  vector<T> B;  // A * Phi
  {
    ifstream value_stream("/home/k/B");
    read_mat(value_stream, K, N, &B);
  }

  vector<T> y;  // measuremrnts
  {
    ifstream value_stream("/home/k/y");
    read_mat(value_stream, K, 1, &y);
  }

  amp(B, y);
  return 0;
}

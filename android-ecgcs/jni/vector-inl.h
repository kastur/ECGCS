#ifndef COMMON_INL_H_
#define COMMON_INL_H_
#include <stdlib.h>
#include <string.h>


template<class T>
class vector {
 public:
  vector(size_t n) : data_(NULL) {
		n_ = n;
		data_ = (T*)calloc(n_, sizeof(T));
		if (data_ == NULL) {
      exit(-1);
    }
	}

  vector(size_t n, const T* data) : data_(NULL) {
		n_ = n;
		data_ = (T*)calloc(n_, sizeof(T));
		if (data_ == NULL) {
      exit(-1);
    }
		memcpy(data_, data, n_*sizeof(T));
	}

	~vector() {
		free(data_);
	}

	void fill(T val) {
		for (int i = 0; i < n_; i++)
			data_[i] = val;
	}

	size_t size() const {
		return n_;
	}

	const T& operator[](const int i) const {
		return data_[i];
	}

	T& operator[](const int i) {
		return data_[i];
	}

	vector<T>& operator=(const vector<T>& rhs) {
		if (this == &rhs)
			return *this;
		memcpy(data_, rhs.data_, rhs.size()*sizeof(T));
		return *this;
	}

 private:
	T* data_;  // array to store elements
	size_t n_;  // size of array
};


#endif  // COMMON_INL_H_

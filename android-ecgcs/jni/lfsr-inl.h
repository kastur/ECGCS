#ifndef LFSR_INL_H_
#define LFSR_INL_H_

#include "vector-inl.h"

class LFSRRand {
 public:
  LFSRRand(unsigned int state = 1, unsigned int coef = 0xd0000001u) :
		state_(state),
		coef_(coef) {};

	unsigned int getState() {
		return state_;
	}

	unsigned int random() {
		 // taps: 32 31 29 1; 
		 // characteristic polynomial: x^32 + x^31 + x^29 + x + 1
		 state_ = (state_ >> 1) ^ (unsigned int)(0 - (state_ & 1u) & coef_); 
		 return state_ & 1u;
	}

 private:
  unsigned int state_;
  unsigned int coef_;
};

template <class T>
unsigned int BernoulliNormalized(LFSRRand& rand, T norm_factor, vector<T>& x) {
	const T one = 1;
  for (int i = 0; i < x.size(); i++) {
    //x[i] = (rand.random() ? one : -one) / norm_factor;
    x[i] = rand.random() / norm_factor;
  }
}

#endif  // LFSR_INL_H_

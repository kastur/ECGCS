#ifndef LARGEST_ELEMENT_INL_H_
#define LARGEST_ELEMENT_INL_H_

#include <stdlib.h>

template<class T>
static int compare(const void* a, const void* b) {
	const T* aval = static_cast<const T*>(a);
	const T* bval = static_cast<const T*>(b);
	if (*aval < *bval) return -1;
	else if (*aval > *bval) return 1;
	else return 0;
}

template <class T>
void largestElement(const vector<T>& signal, size_t k, T& threshold) {
	vector<T> signal_copy(signal.size());
	signal_copy = signal;
	qsort(&signal_copy[0], signal_copy.size(), sizeof(T), compare<T>);
	threshold = signal_copy[signal_copy.size() - k];
}

#endif

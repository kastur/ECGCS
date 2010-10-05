// Author: Kasturi Rangan Raghavan
// \file
// Emulates barebones AnalogIn from mbed library

#ifndef ECG_READER_READER_INL_H_
#define ECG_READER_READER_INL_H_
#include <string>
#include <fstream>

using std::ifstream;
using std::string;

class AnalogIn {
 public:
  AnalogIn(const string& file, float scale, float offset) {
    stream_.open(file.c_str());
    scale_ = scale;
    offset_ = offset;
  }

  float read() {
    float t, mV;
    if(stream_.bad()) {
      return 0.0;
    } else {
      stream_ >> t >> mV;
      return scale_ * mV + offset_;
    }
  }
 private:
  ifstream stream_;
  float scale_, offset_;
};


#endif  // ECG_READER_READER_INL_H_

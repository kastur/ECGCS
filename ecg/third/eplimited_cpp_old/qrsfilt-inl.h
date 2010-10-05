#ifndef ECG_THIRD_EPLIMITED_QRSFILT_INL_H_
#define ECG_THIRD_EPLIMITED_QRSFILT_INL_H_
// Modified to OOP class sturcture
// Author: Kasturi Rangan Raghavan
/*****************************************************************************
FILE:  qrsfilt.cpp
AUTHOR:    Patrick S. Hamilton
REVISED:    5/13/2002
  ___________________________________________________________________________

qrsfilt.cpp filter functions to aid beat detecton in electrocardiograms.
Copywrite (C) 2000 Patrick S. Hamilton

This file is free software; you can redistribute it and/or modify it under
the terms of the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Library General Public License for more
details.

You should have received a copy of the GNU Library General Public License along
with this library; if not, write to the Free Software Foundation, Inc., 59
Temple Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the author by e-mail (pat@eplimited.edu) or postal mail
(Patrick Hamilton, E.P. Limited, 35 Medford St., Suite 204 Somerville,
MA 02143 USA).  For updates to this software, please visit our website
(http://www.eplimited.com).
  __________________________________________________________________________

    This file includes QRSFilt() and associated filtering files used for QRS
    detection.  Only QRSFilt() and deriv1() are called by the QRS detector
    other functions can be hidden.
    Revisions:
        5/13: Filter implementations have been modified to allow simplified
            modification for different sample rates.
*******************************************************************************/
#include <cstdlib>

namespace qrsdet {

/******************************************************************************
* Syntax:
*    int QRSFilter(int datum, int init) ;
* Description:
*    QRSFilter() takes samples of an ECG signal as input and returns a sample of
*    a signal that is an estimate of the local energy in the QRS bandwidth.  In
*    other words, the signal has a lump in it whenever a QRS complex, or QRS
*    complex like artifact occurs.  The filters were originally designed for data
*  sampled at 200 samples per second, but they work nearly as well at sample
*    frequencies from 150 to 250 samples per second.
*
*    The filter buffers and static variables are reset if a value other than
*    0 is passed to QRSFilter through init.
*******************************************************************************/

class RingBuffer {
 public:
  RingBuffer(size_t buf_size) {
    this->buf_size = buf_size;
    data = new int[buf_size];
      for (ptr = 0; ptr < buf_size; ptr++)
          data[ptr] = 0;
    ptr = 0;
  }

  ~RingBuffer() {
    delete[] data;
  }

  int* GetRange(int begin, int end, int* n) {
    *n = end - begin;
    int* mem = new int[*n];
    for (int p = begin, i = 0; p < end; p++, i++) {
      mem[i] = GetValue(p);
    }
    return mem;
  }

  void SetValue(int value, int offset = 0) {
    int set_ptr = ptr + offset;
    if (set_ptr < 0) 
      set_ptr += static_cast<int>(buf_size);
    set_ptr %= buf_size;
    data[set_ptr] = value;
  }

  int GetValue(int offset = 0) {
    int get_ptr = ptr + offset;
    if (get_ptr < 0)
      get_ptr += static_cast<int>(buf_size);
    get_ptr %= buf_size;
    return data[get_ptr];
  }
    
  void Increment() {
    if (++ptr == buf_size)
      ptr = 0;
  }

 protected:
  size_t buf_size;
  int* data;
  int ptr;
};

class DelayFilter : public RingBuffer {
 public:
  explicit DelayFilter(size_t buf_size) : RingBuffer(buf_size) {}
  int Process(int datum) {
    int output = GetValue();
    SetValue(datum);
    Increment();
    return output;
  }
};

/*************************************************************************
*  lpfilt() implements the digital filter represented by the difference
*  equation:
*
*     y[n] = 2*y[n-1] - y[n-2] + x[n] - 2*x[t-24 ms] + x[t-48 ms]
*
*    Note that the filter delay is (LPBUFFER_LGTH/2)-1
*
**************************************************************************/
class LowPassFilter : public RingBuffer {
 public:
  LowPassFilter(size_t buf_size) : RingBuffer(buf_size) {
    y1 = 0;
    y2 = 0;
  }

  int Process(int datum) {
      int length = static_cast<int>(buf_size);
      long y0;
      int output, halfPtr;
      halfPtr = ptr - length / 2;
      if (halfPtr < 0)  // to x[n-6].
      halfPtr += length;
      y0 = (y1 << 1) - y2 + datum - (data[halfPtr] << 1) + data[ptr];
      y2 = y1;
      y1 = y0;
      output = y0 / ((length * length) / 4);
      data[ptr] = datum;
    Increment();
      return output;
  }

 private:
  long y1, y2;
};


/******************************************************************************
*  hpfilt() implements the high pass filter represented by the following
*  difference equation:
*
*    y[n] = y[n-1] + x[n] - x[n-128 ms]
*    z[n] = x[n-64 ms] - y[n] ;
*
*  Filter delay is (buf_size-1)/2
******************************************************************************/
class HighPassFilter : public RingBuffer {
 public:
  HighPassFilter(size_t buf_size) : RingBuffer(buf_size) {
        y = 0;
  }

  int Process(int datum) {
      int output, halfPtr;
      int length = static_cast<int>(buf_size);
    y += datum - data[ptr];
      halfPtr = ptr - (length / 2);
      if (halfPtr < 0)
          halfPtr += length;
      output = data[halfPtr] - (y / length);
      data[ptr] = datum;
      Increment();
      return output;
  }

 private:
  long y;
};


/*****************************************************************************
*  deriv1 and deriv2 implement derivative approximations represented by
*  the difference equation:
*
*    y[n] = x[n] - x[n - 10ms]
*
*  Filter delay is DERIV_LENGTH/2
*****************************************************************************/
class DerivFilter : public RingBuffer {
 public:
  DerivFilter(size_t buf_size) : RingBuffer(buf_size) {}

  int Process(int datum) {
    int y;
    y = datum - data[ptr];
      data[ptr] = datum;
    Increment();
      return y;
  }
};


/*****************************************************************************
* mvwint() implements a moving window integrator.  Actually, mvwint() averages
* the signal values over the last WINDOW_WIDTH samples.
*****************************************************************************/
class MovingWindowAverager : public RingBuffer {
 public:
  MovingWindowAverager(size_t buf_size) : RingBuffer(buf_size) {  
    sum = 0;
  }

  int Process(int datum) {
    int output;
    int length = static_cast<int>(buf_size);
      sum += datum;
      sum -= data[ptr];
      data[ptr] = datum;
    Increment();
      if((sum / length) > 32000)
          output = 32000 ;
      else
          output = sum / length;
      return output;
  }
 private:
  long sum;
};

struct QRSFilterParams {
  QRSFilterParams(double ts) {
    this->ts = ts;
  }

  size_t get_ms(size_t t) const { return static_cast<double>(t) / ts + 0.5; }

  double get_ts() const { return ts; }

  size_t get_lp_length() const { return get_ms(50); }

  size_t get_hp_length() const { return get_ms(125); }

  size_t get_d_length() const { return get_ms(10); }

  size_t get_window_length() const { return get_ms(80); }

  size_t get_preblank() const { return get_ms(200); }

  size_t get_filter_delay() const {
    return get_d_length() / 2.0 +
           get_lp_length() / 2.0 - 1.0 +
           (get_hp_length() - 1.0) / 2.0 +
           get_preblank();
  }

  size_t get_der_delay() const {
    return get_filter_delay() + get_window_length() + get_ms(100);
  }

 private:
  double ts;
};

class QRSFilter {
 public:
  QRSFilter(QRSFilterParams params) :
  lp_filter(params.get_lp_length()),
  hp_filter(params.get_hp_length()),
  d_filter(params.get_d_length()),
  mwa_filter(params.get_window_length()) {}

  int Process(int datum) {
    datum = lp_filter.Process(datum);
    datum = hp_filter.Process(datum);
    datum = d_filter.Process(datum);
    datum = abs(datum);
    datum = mwa_filter.Process(datum);
    return datum;
  }

 private:
  LowPassFilter lp_filter;
  HighPassFilter hp_filter;
  DerivFilter d_filter;
  MovingWindowAverager mwa_filter;
};

/**************************************************************
* peak() takes a datum as input and returns a peak height
* when the signal returns to half its peak height, or 
**************************************************************/
class PeakDetector {
 public:
  PeakDetector(const QRSFilterParams& params) {
    max = 0;
    timeSinceMax = 0;
    lastDatum = 0;
    ms_95_ = params.get_ms(95);
  }

  // what is Dly used for?
  int Process(int datum) {
      int pk = 0;
      if(timeSinceMax > 0)
          ++timeSinceMax ;

      if((datum > lastDatum) && (datum > max)) {
          max = datum;
          if(max > 2)
              timeSinceMax = 1;
        } else if(datum < (max >> 1)) {
      pk = max ;
      max = 0 ;
      timeSinceMax = 0 ;
          Dly = 0 ;
        } else if(timeSinceMax > ms_95_) {
          pk = max;
          max = 0;
          timeSinceMax = 0;
          Dly = 3;
        }
      lastDatum = datum ;
      return pk;
  }
 private:
  size_t ms_95_;
  int Dly;  // unused?
  int max;
  int timeSinceMax;
  int lastDatum;
};

/********************************************************************
mean returns the mean of an array of integers.  It uses a slow
sort algorithm, but these arrays are small, so it hardly matters.
********************************************************************/
int mean(int *array, int n)
    {
    long sum = 0;
    for(int i = 0; i < n; ++i)
        sum += array[i];
    sum /= n;
    return sum;
    }

/****************************************************************************
 thresh() calculates the detection threshold from the qrs mean and noise
 mean estimates.
****************************************************************************/
int detection_thresh(int qmean, int nmean) {
  static const float TH = .3125;
  int thrsh, dmed;
  float temp;
  dmed = qmean - nmean;
  temp = dmed;
  temp *= TH;
  dmed = temp;
  thrsh = nmean + dmed;
  return thrsh;
}

/***********************************************************************
    BLSCheck() reviews data to see if a baseline shift has occurred.
    This is done by looking for both positive and negative slopes of
    roughly the same magnitude in a 220 ms window.
***********************************************************************/
class BLSCheckBuffer : RingBuffer {
 public:
  BLSCheckBuffer(const QRSFilterParams& params)
    : RingBuffer(params.get_der_delay()) {
    ms150_ = params.get_ms(150);
    ms220_ = params.get_ms(220);
  }

  int Process(int datum) {
     data[ptr] = datum;
     Increment();
     return datum;
  }

  int BLSCheck(int *maxder) {
    int max = 0, min = 0;
    int maxt, mint;
    int saved_ptr = ptr;
      for(int t = 0; t < ms220_; t++) {
          
          int x = data[ptr];
          if(x > max) {
              maxt = t;
              max = x;
          } else if(x < min) {
              mint = t;
              min = x;
          }
        Increment();
      }
      ptr = saved_ptr;
      *maxder = max;
      min = -min;
    
      /* Possible beat if a maximum and minimum pair are found
          where the interval between them is less than 150 ms. */
         
      if((max > (min >> 3)) && (min > (max >> 3)) &&
          (abs(maxt - mint) < ms150_))
          return 0;
      else
          return 1;
  }
 private:
  size_t ms150_, ms220_;
};


}  // namespace qrsdet

#endif  // ECG_THIRD_EPLIMITED_QRSFILT_INL_H_


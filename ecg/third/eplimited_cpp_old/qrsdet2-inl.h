#ifndef ECG_THIRD_EPLIMITED_QRSDET2_INL_H_
#define ECG_THIRD_EPLIMITED_QRSDET2_INL_H_
// Modified to OOP class sturcture
// Author: Kasturi Rangan Raghavan
/*****************************************************************************
FILE:  qrsdet2.cpp
AUTHOR:    Patrick S. Hamilton
REVISED:    7/08/2002
  ___________________________________________________________________________

qrsdet2.cpp: A QRS detector.
Copywrite (C) 2002 Patrick S. Hamilton

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

This file contains functions for detecting QRS complexes in an ECG.  The
QRS detector requires filter functions in qrsfilt.cpp and parameter
definitions in qrsdet.h.  QRSDet is the only function that needs to be
visable outside of these files.

Syntax:
    int QRSDet(int ecgSample, int init) ;

Description:
    QRSDet() implements a modified version of the QRS detection
    algorithm described in:

    Hamilton, Tompkins, W. J., "Quantitative investigation of QRS
    detection rules using the MIT/BIH arrhythmia database",
    IEEE Trans. Biomed. Eng., BME-33, pp. 1158-1165, 1987.

    Consecutive ECG samples are passed to QRSDet.  QRSDet was
    designed for a 200 Hz sample rate.  QRSDet contains a number
    of static variables that it uses to adapt to different ECG
    signals.  These variables can be reset by passing any value
    not equal to 0 in init.

    Note: QRSDet() requires filters in QRSFilt.cpp

Returns:
    When a QRS complex is detected QRSDet returns the detection delay.

****************************************************************/

/* For memmove. */
#ifdef __STDC__
#include <string.h>
#else
#include <mem.h>
#endif

#include "qrsfilt-inl.h"


namespace qrsdet {

int mean(int *array, int datnum);

int detection_thresh(int qmean, int nmean);


// TODO(krr) : Where to put this?
const int MEMMOVELEN = 7*sizeof(int);


class QRSDet2 {
 public:
  QRSDet2(const QRSFilterParams& params)
    : qrs_filter(params),
      peak_detector(params),
      first_derivative_filter(params.get_d_length()),
      bls_check_buffer(params) {
    rsetCount = 0;
    sbcount = params.get_ms(1500);
    ms_preblank_ = params.get_ms(200);
    ms_360_ = params.get_ms(360);
    ms_1000_ = params.get_ms(1000);
    ms_1650_ = params.get_ms(1650);
    ms_window_length_ = params.get_window_length();
    ms_filter_delay_ = params.get_filter_delay();

    for (int i = 0; i < 8; i++) {
            noise[i] = 0;    /* Initialize noise buffer */
            rrbuf[i] = ms_1000_;/* and R-to-R interval buffer. */
    }
        qpkcnt = maxder = lastmax = count = sbpeak = 0 ;
        initBlank = initMax = preBlankCnt = 0 ;
  }

  int Process(int datum) {
      int fdatum, QrsDelay = 0;
      int newPeak, aPeak;
      fdatum = qrs_filter.Process(datum);
      /* Wait until normal detector is ready before calling early detections. */
      aPeak = peak_detector.Process(fdatum);
      if (aPeak < QRSDet2::MIN_PEAK_AMP)
          aPeak = 0;

      // Hold any peak that is detected for 200 ms
      // in case a bigger one comes along. There
      // can only be one QRS complex in any 200 ms window.

      newPeak = 0;
      if (aPeak && !preBlankCnt) {
      // If there has been no peak for 200 ms
      // save this one and start counting.
          tempPeak = aPeak;
          preBlankCnt = ms_preblank_;
        } else if (!aPeak && preBlankCnt) {          
      // If we have held onto a peak for
      // 200 ms pass it on for evaluation.
          if (--preBlankCnt == 0)
              newPeak = tempPeak ;
        } else if (aPeak) {
      // If we were holding a peak, but
      // this ones bigger, save it and
          if (aPeak > tempPeak) {
        // start counting to 200 ms again.
              tempPeak = aPeak ;
        preBlankCnt = ms_preblank_;
      } else if (--preBlankCnt == 0) {
              newPeak = tempPeak ;
      }
        }

      // Save derivative of raw signal for T-wave and baseline
      //   shift discrimination.
      int deriv_datum = first_derivative_filter.Process(datum);
    bls_check_buffer.Process(deriv_datum);

      // Initialize the qrs peak buffer with the first eight
      // local maximum peaks detected.
      if (qpkcnt < 8) {
          ++count;
          if (newPeak > 0)
        count = ms_window_length_;
          if (++initBlank == ms_1000_) {
              initBlank = 0;
              qrsbuf[qpkcnt] = initMax;
              initMax = 0;
              ++qpkcnt;
              if (qpkcnt == 8) {
                  qmean = mean(qrsbuf, 8) ;
                  nmean = 0 ;
                  rrmean = ms_1000_ ;
                  sbcount = ms_1650_ ;
                  det_thresh = detection_thresh(qmean,nmean) ;
                  }
              }
          if (newPeak > initMax)
              initMax = newPeak ;
          }

      else    /* Else test for a qrs. */
          {
          ++count ;
          if(newPeak > 0)
              {
            
            
              /* Check for maximum derivative and matching minima and maxima
                 for T-wave and baseline shift rejection.  Only consider this
                 peak if it doesn't seem to be a base line shift. */
                 
              if(!bls_check_buffer.BLSCheck(&maxder))
                  {


                  // Classify the beat as a QRS complex
                  // if the peak is larger than the detection threshold.

                  if(newPeak > det_thresh)
                      {
                      memmove(&qrsbuf[1], qrsbuf, MEMMOVELEN) ;
                      qrsbuf[0] = newPeak ;
                      qmean = mean(qrsbuf,8) ;
                      det_thresh = detection_thresh(qmean,nmean) ;
                      memmove(&rrbuf[1], rrbuf, MEMMOVELEN) ;
                      rrbuf[0] = count - ms_window_length_ ;
                      rrmean = mean(rrbuf,8) ;
                      sbcount = rrmean + (rrmean >> 1) + ms_window_length_ ;
                      count = ms_window_length_ ;

                      sbpeak = 0 ;

                      lastmax = maxder ;
                      maxder = 0 ;
                      QrsDelay =  ms_window_length_ + ms_filter_delay_ ;
                      initBlank = initMax = rsetCount = 0 ;
                      }

                  // If a peak isn't a QRS update noise buffer and estimate.
                  // Store the peak for possible search back.


                  else
                      {
                      memmove(&noise[1],noise,MEMMOVELEN) ;
                      noise[0] = newPeak ;
                      nmean = mean(noise,8) ;
                      det_thresh = detection_thresh(qmean,nmean) ;

                      // Don't include early peaks (which might be T-waves)
                      // in the search back process.  A T-wave can mask
                      // a small following QRS.

                      if((newPeak > sbpeak) && ((count-ms_window_length_) >= ms_360_))
                          {
                          sbpeak = newPeak ;
                          sbloc = count  - ms_window_length_ ;
                          }
                      }
                  }
              }
        
          /* Test for search back condition.  If a QRS is found in  */
          /* search back update the QRS buffer and det_thresh.      */

          if((count > sbcount) && (sbpeak > (det_thresh >> 1)))
              {
              memmove(&qrsbuf[1],qrsbuf,MEMMOVELEN) ;
              qrsbuf[0] = sbpeak ;
              qmean = mean(qrsbuf,8) ;
              det_thresh = detection_thresh(qmean,nmean) ;
              memmove(&rrbuf[1],rrbuf,MEMMOVELEN) ;
              rrbuf[0] = sbloc ;
              rrmean = mean(rrbuf,8) ;
              sbcount = rrmean + (rrmean >> 1) + ms_window_length_ ;
              QrsDelay = count = count - sbloc ;
              QrsDelay += ms_filter_delay_;
              sbpeak = 0 ;
              lastmax = maxder ;
              maxder = 0 ;

              initBlank = initMax = rsetCount = 0 ;
              }
          }

      // In the background estimate threshold to replace adaptive threshold
      // if eight seconds elapses without a QRS detection.

      if( qpkcnt == 8 )
          {
          if(++initBlank == ms_1000_)
              {
              initBlank = 0 ;
              rsetBuff[rsetCount] = initMax ;
              initMax = 0 ;
              ++rsetCount ;

              // Reset threshold if it has been 8 seconds without
              // a detection.

              if(rsetCount == 8)
                  {
                  for(int i = 0; i < 8; ++i)
                      {
                      qrsbuf[i] = rsetBuff[i] ;
                      noise[i] = 0 ;
                      }
                  qmean = mean( rsetBuff, 8 ) ;
                  nmean = 0 ;
                  rrmean = ms_1000_;
                  sbcount = ms_1650_;
                  det_thresh = detection_thresh(qmean,nmean) ;
                  initBlank = initMax = rsetCount = 0 ;
                  }
              }
          if( newPeak > initMax )
              initMax = newPeak ;
          }

      return(QrsDelay);
  }
 private:
  int det_thresh, qpkcnt;
    int qrsbuf[8], noise[8], rrbuf[8];
    int rsetBuff[8], rsetCount;
    int nmean, qmean, rrmean;
    int count, sbpeak, sbloc, sbcount;
    int maxder, lastmax;
    int initBlank, initMax;
    int preBlankCnt, tempPeak;
  QRSFilter qrs_filter;
  PeakDetector peak_detector;
  DerivFilter first_derivative_filter;
  BLSCheckBuffer bls_check_buffer;
  size_t ms_preblank_;
  size_t ms_360_;
  size_t ms_1000_;
  size_t ms_1650_;
  size_t ms_window_length_;
  size_t ms_filter_delay_;

  // Prevents detections of peaks smaller than 150 uV.
  static const int MIN_PEAK_AMP =  7;
};


}  // namespace qrsdet

#endif  // ECG_THIRD_EPLIMITED_QRSDET2_INL_H_



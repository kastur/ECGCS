#ifndef ECG_THIRD_EPLIMITED_QRSDET_H_
#define ECG_THIRD_EPLIMITED_QRSDET_H_

/*****************************************************************************
FILE:  qrsdet.h
AUTHOR:	Patrick S. Hamilton
REVISED:	4/16/2002
  ___________________________________________________________________________

qrsdet.h QRS detector parameter definitions
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

You may contact the author by e-mail (pat@eplimited.com) or postal mail
(Patrick Hamilton, E.P. Limited, 35 Medford St., Suite 204 Somerville,
MA 02143 USA).  For updates to this software, please visit our website
(http://www.eplimited.com).
  __________________________________________________________________________
  Revisions:
	4/16: Modified to allow simplified modification of digital filters in
   	qrsfilt().
*****************************************************************************/
/*
namespace qrsdet {

static const int SAMPLE_RATE	= 125;  // Sample rate in Hz.
static const double MS_PER_SAMPLE = 1000.0 / static_cast<double>(SAMPLE_RATE);

static const int MS10 = 10.0/MS_PER_SAMPLE + 0.5;
static const int MS25 = 25.0/MS_PER_SAMPLE + 0.5;
static const int MS30 = 30.0/MS_PER_SAMPLE + 0.5;
static const int MS80 = 80.0/MS_PER_SAMPLE + 0.5;
static const int MS95 = 95.0/MS_PER_SAMPLE + 0.5;
static const int MS100 = 100.0/MS_PER_SAMPLE + 0.5;
static const int MS125	= 125.0/MS_PER_SAMPLE + 0.5;
static const int MS150	= 150.0/MS_PER_SAMPLE + 0.5;
static const int MS160	= 160.0/MS_PER_SAMPLE + 0.5;
static const int MS175	= 175.0/MS_PER_SAMPLE + 0.5;
static const int MS195	= 195.0/MS_PER_SAMPLE + 0.5;
static const int MS200	= 200.0/MS_PER_SAMPLE + 0.5;
static const int MS220	= 220.0/MS_PER_SAMPLE + 0.5;
static const int MS250	= 250.0/MS_PER_SAMPLE + 0.5;
static const int MS300	= 300.0/MS_PER_SAMPLE + 0.5;
static const int MS360	= 360.0/MS_PER_SAMPLE + 0.5;
static const int MS450	= 450.0/MS_PER_SAMPLE + 0.5;
static const int MS1000 = SAMPLE_RATE;
static const int MS1500 = 1500.0/MS_PER_SAMPLE;

static const int PRE_BLANK = MS200;

// Filter delays plus 200 ms blanking delay

}  // end namespace qrsdet
*/
#endif  // ECG_THIRD_EPLIMITED_QRSDET_H_


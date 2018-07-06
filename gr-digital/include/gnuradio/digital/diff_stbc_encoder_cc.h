/* -*- c++ -*- */
/* 
 * Copyright 2018 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifndef INCLUDED_DIGITAL_DIFF_STBC_ENCODER_CC_H
#define INCLUDED_DIGITAL_DIFF_STBC_ENCODER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace digital {

/*! \brief Encodes an incoming stream after the rules of a differential
 * Space Time Block Code (STBC), producing two output streams.
 *
 * The differential STBC encoder works with sequences of length 2, encoding them
 * into a code sequence of length 2 at each of the two output branches, respectively.
 * The number of input items is automatically scheduled to a multiple of 2, however,
 * if the input data stream is terminated, the absolute number of input items
 * must be an even number.
 * The differential STBC encoder is a sync block which produces the same amount of
 * output items at each port as there are input items. The code rate is R=1.
 *
 * The incoming samples are expected to be PSK modulated samples of any modulation order.
 * If the constellation is phase shifted, meaning that there is no constellation point
 * on the real axis, this must be stated by setting the input argument 'phase_offset'
 * to the phase (in radians) of one of the constellation points.
 *
 * The algorithm of this differential STBC follows [1]. Briefly, this STBC
 * is the differential version of the well known Alamouti Code.
 *
 * [1] V. Tarokh and H. Jafarkhani, "A differential detection scheme for transmit diversity,"
 * WCNC. 1999 IEEE Wireless Communications and Networking Conference (Cat. No.99TH8466),
 * New Orleans, LA, 1999, pp. 1043-1047 vol.3. doi: 10.1109/WCNC.1999.796832
 */
    class DIGITAL_API diff_stbc_encoder_cc : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<diff_stbc_encoder_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::diff_stbc_encoder_cc.
       *
       * To avoid accidental use of raw pointers, digital::diff_stbc_encoder_cc's
       * constructor is in a private implementation
       * class. digital::diff_stbc_encoder_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(float phase_offset = 0.0);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIFF_STBC_ENCODER_CC_H */


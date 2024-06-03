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


#ifndef INCLUDED_DIGITAL_ALAMOUTI_ENCODER_CC_H
#define INCLUDED_DIGITAL_ALAMOUTI_ENCODER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace digital {

/*! \brief Encodes an incoming stream after the rules of Alamouti's code,
 * producing two output streams.
 *  \ingroup digital
 *
 * The Alamouti encoder works with sequences of length 2, encoding them
 * into a code sequence of length 2 at each of the two output branches, respectively.
 * The number of input items is automatically scheduled to a multiple of 2, however,
 * if the input data stream is terminated, the absolute number of input items
 * must be an even number.
 * The Alamouti encoder is a sync block which produces the same amount of
 * output items at each port as there are input items. The code rate is R=1.
 *
 * There exist different versions of the exact algorithm (changed
 * position of negations). This implementation follows [1].
 *
 * [1] Andrea Goldsmith. 2005. Wireless Communications.
 *     Cambridge University Press, New York, NY, USA.
 */
    class DIGITAL_API alamouti_encoder_cc : virtual public gr::sync_block
    {
     public:
      typedef std::shared_ptr<alamouti_encoder_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::alamouti_encoder_cc.
       *
       * To avoid accidental use of raw pointers, digital::alamouti_encoder_cc's
       * constructor is in a private implementation
       * class. digital::alamouti_encoder_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint32_t vlen=1);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_ALAMOUTI_ENCODER_CC_H */


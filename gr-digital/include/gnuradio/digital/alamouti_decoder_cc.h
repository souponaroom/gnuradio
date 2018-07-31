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


#ifndef INCLUDED_DIGITAL_ALAMOUTI_DECODER_CC_H
#define INCLUDED_DIGITAL_ALAMOUTI_DECODER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_interpolator.h>

namespace gr {
  namespace digital {

/*! \brief Decodes an incoming stream after the rules of Alamouti's code.
 *  \ingroup digital
 *
 * The Alamouti decoder works with sequences of length 2, decoding them
 * into a code sequence of length 2.
 * The number of input items is automatically scheduled to a multiple of 2, however,
 * if the input data stream is terminated, the absolute number of input items
 * must be an even number.
 * The Alamouti decoder is a sync block which produces the same amount of
 * output items as there are input items. The code rate is R=1.
 *
 * The CSI is transported via stream tags with key='csi'.
 * Initially the CSI is set to 1.0 + 0j for both branches and are updated
 * with each incoming CSI. Because the Alamouti algorithm works with sequences
 * of length 2, the tags should be set only on samples with even positions.
 * CSI tags on uneven sample positions are not processed until the beginning of
 * the next sequence (in this case one sample delay) begins.
 *
 * There exist different versions of the exact algorithm (changed
 * position of negations). This implementation follows [1] and is therefore
 * consistent with the Alamouti encoding block 'alamouti_encoder_cc'.
 *
 * [1] Andrea Goldsmith. 2005. Wireless Communications.
 *     Cambridge University Press, New York, NY, USA.
 */
    class DIGITAL_API alamouti_decoder_cc : virtual public gr::sync_interpolator
    {
     public:
      typedef boost::shared_ptr<alamouti_decoder_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::alamouti_decoder_cc.
       *
       * To avoid accidental use of raw pointers, digital::alamouti_decoder_cc's
       * constructor is in a private implementation
       * class. digital::alamouti_decoder_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint32_t vlen=1);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_ALAMOUTI_DECODER_CC_H */


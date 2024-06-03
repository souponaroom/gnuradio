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


#ifndef INCLUDED_DIGITAL_DIFF_STBC_DECODER_CC_H
#define INCLUDED_DIGITAL_DIFF_STBC_DECODER_CC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace digital {

/*! \brief Decodes an incoming stream after the rules of a differential
 * Space Time Block Code (STBC).
 *
 * The differential STBC decoder works with 'sequences' of length 2, decoding them
 * into a code sequence of length 2.
 * The number of input items is automatically scheduled to a multiple of 2, however,
 * if the input data stream is terminated, the absolute number of input items
 * must be an even number.
 * When running with an infinite sequence stream, the differential STBC decoder is a sync
 * block which produces the same amount of output items as there are input items.
 * The code rate is R=1.
 * If a finite block of sequences is transmitted, the decoder needs a reference sequence
 * to differentially decode the first sequence. When using stream tags 'start'
 * on the stream, the decoder therefore produces one sequence less than it consumes, because
 * the tagged sequence is skipped and taken as reference to differentially decode the first
 * sequence after the tag as first decoded output sequence. If there is no tag at the very
 * first sequence of the transmission, the decoder uses a dummy sequence (= basis vectors of
 * basis transformation 'd_basis_vecs') to decode the first sequence. Generally this leads
 * to a wrong decoded first symbol. The following symbols are decoded correctly.
 * To avoid a wrong decoded first symbol, there must be a stream tag on the first
 * incoming sequence.
 * You can set a correct stream tag in a C++ block like in the following example,
 * where offset stands for the item position of the tag in the current output buffer of a call of the work function:
 * add_item_tag(0, nitems_written(0) + offset, pmt::mp("start"), 0);
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
    class DIGITAL_API diff_stbc_decoder_cc : virtual public gr::block
    {
     public:
      typedef std::shared_ptr<diff_stbc_decoder_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::diff_stbc_decoder_cc.
       *
       * To avoid accidental use of raw pointers, digital::diff_stbc_decoder_cc's
       * constructor is in a private implementation
       * class. digital::diff_stbc_decoder_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(float phase_offset=0.0, uint32_t vlen=1, const std::string &start_key="start");
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIFF_STBC_DECODER_CC_H */


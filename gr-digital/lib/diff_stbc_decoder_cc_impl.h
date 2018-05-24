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

#ifndef INCLUDED_DIGITAL_DIFF_STBC_DECODER_CC_IMPL_H
#define INCLUDED_DIGITAL_DIFF_STBC_DECODER_CC_IMPL_H

#include <gnuradio/digital/diff_stbc_decoder_cc.h>

namespace gr {
  namespace digital {

    class diff_stbc_decoder_cc_impl : public diff_stbc_decoder_cc
    {
     private:
      std::vector <gr::tag_t> tags; /*!< Vector that stores the tags in input buffer. */
      static const std::string s; /*!< String that matches the key of the CSI tags. */
      static const pmt::pmt_t d_key; /*!< PMT stores the key of the CSI tag. */
      const std::vector <gr_complex> d_basis_vecs;
      /*!< Complex vector of size 2. The 2 complex elements of this vectors
       * can each be interpreted as a basis vector in the complex plane.
       * These 2 2-dimensional vectors (= 2 complex numbers) are a new basis
       * for the complex plane in which the input samples at the encoder were
       * transformed into and from which the input samples of the decoder are
       * transformed back to the standard basis {(1,0),(0,1)}.
       * [1] refers to this vector as '(a_1, a_2)'.
       */
      gr_complex d_predecessor[2];
      /*!< Complex array of size 2 which stores the last sequence of the
       * previous input buffer.
       */

      void decode_sequences(const gr_complex* prev_seq,
                            const gr_complex* seq,
                            gr_complex* out,
                            uint32_t length);

     public:
      diff_stbc_decoder_cc_impl(float phase_offset);
      ~diff_stbc_decoder_cc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIFF_STBC_DECODER_CC_IMPL_H */


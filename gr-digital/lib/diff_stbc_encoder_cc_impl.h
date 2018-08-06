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

#ifndef INCLUDED_DIGITAL_DIFF_STBC_ENCODER_CC_IMPL_H
#define INCLUDED_DIGITAL_DIFF_STBC_ENCODER_CC_IMPL_H

#include <gnuradio/digital/diff_stbc_encoder_cc.h>

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
    class diff_stbc_encoder_cc_impl : public diff_stbc_encoder_cc
    {
     private:
      /*!< Number of items which are processed parallel. This equals a
       * vector size, but with stream input/output item sizes.
       */
      uint32_t d_block_len;
      /*!< Complex vector of size 2. The 2 complex elements of this vectors
       * can each be interpreted as a basis vector in the complex plane.
       * These 2 2-dimensional vectors (= 2 complex numbers) are a new basis
       * for the complex plane in which the input samples are transformed into.
       * [1] refers to this vector as '(a_1, a_2)'.
       */
      const std::vector <gr_complex> d_basis_vecs;
      /*!< Complex vector of size 2. The elements are the coefficients which
       * describe the input sample with the new basis vectors d_basis_vecs.
       * [1] refers to this vector as 'M(S) = (A(S), B(S)). Note that [1] includes
       * bit mapping into this mapping whereas this block accepts complex symbols.
       */
      std::vector <std::vector<gr_complex> > d_mapping_coeffs;
      /*!< 2-dimensional, complex array of size 2 x block_len which stores the first element of the last
       * transmitted sequence of the previous buffer.
       */
      gr_complex *d_predecessor;

      /*!< \brief Calculates the coefficients to a new vector basis of a incoming vector.
       * The calculation of this function equals a basis transformation. The incoming
       * vector is not changed. The result is written to the vector d_mapping_coeffs.
       *
       * @param in Complex pointer to input samples.
       */
      void map_symbols(const gr_complex* in);

      /*!< \brief Calculates the output data which is being transmitted.
       * Briefly, the phase difference to the previous output data is calculated by a complex multiplication.
       * The result is distributed to 2 antennas and 2 time slots after the rules of Alamouti.
       *
       * @param in Complex pointer to the input samples.
       * @param predecessor1 Predecessor at output (=antenna) 1. This points to the first sample of the sequence.
       * @param predecessor2 Predecessor at output (=antenna) 2. This points to the first sample of the sequence.
       * @param out1 Complex pointer to first unwritten sample of output port 1.
       * @param out2 Complex pointer to first unwritten sample of output port 2.
       */
      void encode_data(const gr_complex* in,
                       const gr_complex* predecessor1,
                       const gr_complex* predecessor2,
                       gr_complex* out1,
                       gr_complex* out2);

      void encode_data(const gr_complex* in,
                       gr_complex* out1,
                       gr_complex* out2,
                       uint32_t length);

     public:
      diff_stbc_encoder_cc_impl(float phase_offset, uint32_t block_len);
      ~diff_stbc_encoder_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_DIFF_STBC_ENCODER_CC_IMPL_H */


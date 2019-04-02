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

#ifndef INCLUDED_DIGITAL_VBLAST_DECODER_CC_IMPL_H
#define INCLUDED_DIGITAL_VBLAST_DECODER_CC_IMPL_H

#include <gnuradio/digital/vblast_decoder_cc.h>

namespace gr {
  namespace digital {
/*! \brief VBLAST decoder.
 * Decodes incoming NxN MIMO data by equalizing the received signals in order
 * to extract the M different transmission signals. For the equalization,
 * you can choose between a zero forcing (key='ZF') and
 * a minimum mean squared error (key='MMSE') scheme. Both schemes require
 * channel state information and MMSE additionally needs SNR information.
 *
 * The CSI and SNR information is transported via stream tags with keys given
 * by input arguments.
 * The CSI must be a 3-dimensional vector with the dimensions vlen, num_inputs, num_inputs.
 * For vlen>1 there is provided a separate channel matrix for each
 * vector element. Set vlen>1 if you use a multicarrier system like OFDM
 * (in this case it would be vlen=number of occupied sub-carriers).
 * To generate a proper CSI tag, use the pmt structure pmt_vector(pmt_vector(pmt_c32vector))).
 * Initially the CSI is set to 1.0 + 0j for all branches and are updated
 * with each incoming CSI. The SNR is initially set to 1.0e6; in this case
 * the MMSE equalizer is equal to the ZF equalizer.
 * The CSI and SNR tags are processed separately and they can therefore arrive
 * at different positions of the stream and occur in different frequencies.
 *
 * For 1x1 and 2x2 MIMO schemes, the equalizer is calculated internally.
 * For MxM schemes with M > 2, the C++ template library Eigen is used
 * for the linear algebra operations.
 */
    class vblast_decoder_cc_impl : public vblast_decoder_cc
    {
     private:
      uint16_t d_num_inputs; /*!< Number of inputs ports. This equals the interpolation rate.*/
      /*! Equalization technique zero forcing 'ZF' or minimum mean squared error 'MMSE'. */
      std::string d_equalizer_type;
      uint32_t d_vlen; /*!< Vector length of the incoming items. */
      std::vector <gr::tag_t> tags; /*!< Vector that stores the tags in input buffer. */
      const pmt::pmt_t d_csi_key; /*!< PMT stores the key of the CSI tag. */
      const pmt::pmt_t d_snr_key; /*!< PMT stores the key of the SNR tag. */
      const pmt::pmt_t d_start_key; /*!< PMT stores the key of the <start of frame> tag. */
      std::vector<std::vector<std::vector<gr_complex> > > d_csi; /*!< Current channel matrix. */
      /*! Equalizer matrix.
       * The left-sided matrix multiplication of the mimo_equalizer with the reception vector produces
       * the decoded data vector. */
      std::vector<std::vector<std::vector<gr_complex> > > d_mimo_equalizer;
      std::vector<float> d_snr; /*!< Vector that stores the current SNR values of the M receivers.*/

      /*!< \brief Equalizes the input data with the d_mimo_equalizer matrix. */
      void update_mimo_equalizer();

      /*!< \brief Calculates the new d_mimo_equalizer matrix with the updated CSI or SNR info.
       * The equalization schemes zero forcing ('ZF') or
       * minimum mean squared error ('MMSE') can be chosen.
       * For 1x1 and 2x2 MIMO schemes, the equalizer is calculated internally.
       * For MxM schemes with M > 2, the C++ template library Eigen
       * is used for the linear algebra operations.
       *
       * @param input Vector of pointers to the input items, one entry per input stream.
       * @param out Pointer to the start of unwritten output items.
       * @param offset Number of already written input items. The number of already read items per stream is offset/d_num_inputs.
       * @param length Number of items of the current symbol in the current buffer.
       */
      void equalize_symbol(gr_vector_const_void_star input,
                           gr_complex* out, uint32_t offset,
                           uint32_t length);

     public:
      vblast_decoder_cc_impl(uint16_t num_inputs, std::string equalizer_type,
                             uint32_t vlen, const std::string &csi_tag_key,
                             const std::string &snr_tag_key,
                             const std::string &frame_tag_key);
      ~vblast_decoder_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_VBLAST_DECODER_CC_IMPL_H */


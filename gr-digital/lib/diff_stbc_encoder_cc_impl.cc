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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <boost/format.hpp>
#include <gnuradio/io_signature.h>
#include "diff_stbc_encoder_cc_impl.h"

using namespace boost;

namespace gr {
  namespace digital {

    diff_stbc_encoder_cc::sptr
    diff_stbc_encoder_cc::make(float phase_offset, uint32_t block_len)
    {
      return gnuradio::get_initial_sptr
        (new diff_stbc_encoder_cc_impl(phase_offset, block_len));
    }

    /*
     * The private constructor
     */
    diff_stbc_encoder_cc_impl::diff_stbc_encoder_cc_impl(float phase_offset, uint32_t block_len)
      : gr::sync_block("diff_stbc_encoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(2, 2, sizeof(gr_complex))),
        d_block_len(block_len),
        d_basis_vecs(std::vector<gr_complex>(2, std::polar((float)M_SQRT1_2, phase_offset)))
    {
      d_mapping_coeffs = std::vector<std::vector<gr_complex> > (2, std::vector<gr_complex>(block_len, 0.0));
      d_predecessor = new gr_complex[2*block_len];
      for (unsigned int i = 0; i < 2*block_len; ++i) {
        d_predecessor[i] = d_basis_vecs[0];
      }
      /* Set the number of input and output items to a multiple of 2 blocks,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2*block_len);
    }

    /*
     * Our virtual destructor.
     */
    diff_stbc_encoder_cc_impl::~diff_stbc_encoder_cc_impl()
    {
    }

    void
    diff_stbc_encoder_cc_impl::map_symbols(const gr_complex *in) {
      for (unsigned int k = 0; k < d_block_len; ++k) {
        d_mapping_coeffs[0][k] = in[0+k] * std::conj(d_basis_vecs[0]) + in[1*d_block_len+k] * std::conj(d_basis_vecs[1]);
        d_mapping_coeffs[1][k] = -in[0+k] * d_basis_vecs[1] + in[1*d_block_len+k] * d_basis_vecs[0];
      }
    }

    void
    diff_stbc_encoder_cc_impl::encode_data(const gr_complex *in,
                                           const gr_complex *predecessor1,
                                           const gr_complex *predecessor2,
                                           gr_complex *out1,
                                           gr_complex *out2) {
      // Transform input vector to new basis and calculate new coefficients.
      map_symbols(in);
      for (unsigned int k = 0; k < d_block_len; ++k) { // Iterate over elements of one block.
        // Calculate the output of antenna 1.
        out1[k] = d_mapping_coeffs[0][k] * predecessor1[k] - d_mapping_coeffs[1][k] * std::conj(predecessor2[k]);
        // Calculate the output of antenna 2.
        out2[k] = d_mapping_coeffs[0][k] * predecessor2[k] + d_mapping_coeffs[1][k] * std::conj(predecessor1[k]);
        // Calculate the second element of the output sequence after the rules of Alamouti.
        out1[1*d_block_len+k] = -std::conj(out2[k]);
        out2[1*d_block_len+k] = std::conj(out1[k]);
      }
    }

    void
    diff_stbc_encoder_cc_impl::encode_data(const gr_complex *in,
                                           gr_complex *out1,
                                           gr_complex *out2,
                                           uint32_t length) {
      uint32_t count = 0;

      while (count < length*d_block_len) {
        // Transform input vector to new basis and calculate new coefficients.
        map_symbols(&in[count]);
        for (unsigned int k = 0; k < d_block_len; ++k) {
          // Calculate the output of antenna 1.
          out1[count+2*d_block_len+k] = d_mapping_coeffs[0][k] * out1[count+k] - d_mapping_coeffs[1][k] * std::conj(out2[count+k]);
          // Calculate the output of antenna 2.
          out2[count+2*d_block_len+k] = d_mapping_coeffs[0][k] * out2[count+k] + d_mapping_coeffs[1][k] * std::conj(out1[count+k]);
          // Calculate the second element of the output sequence after the rules of Alamouti.
          out1[count+3*d_block_len+k] = -std::conj(out2[count+2*d_block_len+k]);
          out2[count+3*d_block_len+k] = std::conj(out1[count+2*d_block_len+k]);
        }

        count += 2*d_block_len;
      }
    }

    int
    diff_stbc_encoder_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out1 = (gr_complex *) output_items[0];
      gr_complex *out2 = (gr_complex *) output_items[1];

      // Handle first sequence manually, because of a special predecessor assignment.
      encode_data(in, &d_predecessor[0], &d_predecessor[d_block_len], out1, out2);

      // Calculate the output of both antennas.
      encode_data(&in[2*d_block_len], out1, out2, noutput_items-2*d_block_len);

      // Update predecessor for next call of work.
      if(noutput_items >= 2*d_block_len) {
        for (unsigned int k = 0; k < d_block_len; ++k) {
          d_predecessor[k] = out1[noutput_items - 2*d_block_len + k];
          d_predecessor[d_block_len+k] = out2[noutput_items - 2*d_block_len + k];
        }
      }
      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


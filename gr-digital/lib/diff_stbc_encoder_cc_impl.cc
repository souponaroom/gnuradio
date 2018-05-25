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

#include <gnuradio/io_signature.h>
#include "diff_stbc_encoder_cc_impl.h"

namespace gr {
  namespace digital {

    diff_stbc_encoder_cc::sptr
    diff_stbc_encoder_cc::make(float phase_offset)
    {
      return gnuradio::get_initial_sptr
        (new diff_stbc_encoder_cc_impl(phase_offset));
    }

    /*
     * The private constructor
     */
    diff_stbc_encoder_cc_impl::diff_stbc_encoder_cc_impl(float phase_offset)
      : gr::sync_block("diff_stbc_encoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(2, 2, sizeof(gr_complex))),
        d_basis_vecs(std::vector<gr_complex>(2, std::polar((float)M_SQRT1_2, phase_offset)))
    {
      d_mapping_coeffs = std::vector<gr_complex>(2, 0.0);
      d_predecessor[0] = d_predecessor[1] = d_basis_vecs[0];
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2);
    }

    /*
     * Our virtual destructor.
     */
    diff_stbc_encoder_cc_impl::~diff_stbc_encoder_cc_impl()
    {
    }

    void
    diff_stbc_encoder_cc_impl::map_symbols(const gr_complex *in) {
      d_mapping_coeffs[0] =  in[0]*std::conj(d_basis_vecs[0]) + in[1]*std::conj(d_basis_vecs[1]);
      d_mapping_coeffs[1] = -in[0]*          d_basis_vecs[1]  + in[1]*          d_basis_vecs[0];
    }

    void
    diff_stbc_encoder_cc_impl::calculate_output(const gr_complex *in,
                                        const gr_complex predecessor1,
                                        const gr_complex predecessor2,
                                        gr_complex *out1,
                                        gr_complex *out2) {
      // Transform input vector to new basis and calculate new coefficients.
      map_symbols(in);
      // Calculate the output of antenna 1.
      out1[0] = d_mapping_coeffs[0]*predecessor1 - d_mapping_coeffs[1]*std::conj(predecessor2);
      // Calculate the output of antenna 2.
      out2[0] = d_mapping_coeffs[0]*predecessor2 + d_mapping_coeffs[1]*std::conj(predecessor1);
      // Calculate the second element of the output sequence after the rules of Alamouti.
      out1[1] = -std::conj(out2[0]);
      out2[1] =  std::conj(out1[0]);
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
      calculate_output(in, d_predecessor[0], d_predecessor[1], out1, out2);

      // Iterate over the remaining input sequences.
      for (int i = 1; i < noutput_items/2; ++i) {
        // Calculate the output of both antennas.
        calculate_output(&in[i*2], out1[(i-1)*2], out2[(i-1)*2], &out1[i*2], &out2[i*2]);
      }
      // Update predecessor for next call of work.
      if(noutput_items > 1) {
        d_predecessor[0] = out1[noutput_items - 2];
        d_predecessor[1] = out2[noutput_items - 2];
      }
      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


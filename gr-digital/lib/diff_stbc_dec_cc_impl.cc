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
#include <cmath>
#include <complex>

#include <gnuradio/io_signature.h>
#include "diff_stbc_dec_cc_impl.h"

namespace gr {
  namespace digital {

    diff_stbc_dec_cc::sptr
    diff_stbc_dec_cc::make(float phase_offset)
    {
      return gnuradio::get_initial_sptr
        (new diff_stbc_dec_cc_impl(phase_offset));
    }

    /*
     * The private constructor
     */
    diff_stbc_dec_cc_impl::diff_stbc_dec_cc_impl(float phase_offset)
      : gr::sync_block("diff_stbc_dec_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_basis_vecs(std::vector<gr_complex>(2, std::polar((float)M_SQRT1_2, phase_offset)))
    {
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2);
    }

    /*
     * Our virtual destructor.
     */
    diff_stbc_dec_cc_impl::~diff_stbc_dec_cc_impl()
    {
    }

    int
    diff_stbc_dec_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      // Iterate over the input sequences.
      for (int i = 0; i < (noutput_items/2)-1; ++i) {
        // Decode the received sequence with help of its successor.
        // Calculate the dot product of the received input sequences.
        gr_complex r_1 = in[2*i+2]*std::conj(in[2*i]) + std::conj(in[2*i+3])*in[2*i+1];
        gr_complex r_2 = in[2*i+2]*std::conj(in[2*i+1]) - std::conj(in[2*i+3])*in[2*i];
        // Calculate the decoded (but not normalized) samples and write them to the output buffer.
        out[2*i]   = d_basis_vecs[0] * r_1 - std::conj(d_basis_vecs[1]) * r_2;
        out[2*i+1] = d_basis_vecs[1] * r_1 + std::conj(d_basis_vecs[0]) * r_2;
      }

      // We produced all input sequences except the last one, because it has no successor.
      return noutput_items-2;
    }

  } /* namespace digital */
} /* namespace gr */


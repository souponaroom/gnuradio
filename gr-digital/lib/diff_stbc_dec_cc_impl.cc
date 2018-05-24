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

    const std::string diff_stbc_dec_cc_impl::s = "start";
    const pmt::pmt_t diff_stbc_dec_cc_impl::d_key = pmt::string_to_symbol(s);

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
      // Init predecessor with dummy sequence.
      d_predecessor[0] = d_predecessor[1] = d_basis_vecs[0];
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

    void
    diff_stbc_dec_cc_impl::calculate_output(const gr_complex *in,
                                            const gr_complex predecessor1,
                                            const gr_complex predecessor2,
                                            gr_complex *out) {
      // Calculate the dot product of the received input sequences.
      gr_complex r_1 = in[0]*std::conj(predecessor1) + std::conj(in[1])*predecessor2;
      gr_complex r_2 = in[0]*std::conj(predecessor2) - std::conj(in[1])*predecessor1;
      // Calculate the decoded (but not normalized) samples and write them to the output buffer.
      out[0] = d_basis_vecs[0] * r_1 - std::conj(d_basis_vecs[1]) * r_2;
      out[1] = d_basis_vecs[1] * r_1 + std::conj(d_basis_vecs[0]) * r_2;
    }

    int
    diff_stbc_dec_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      // Collect all tags of the input buffer with key "start" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);
      //TODO check for tags at uneven positions and fix them

      // Iterate over tags.
      for (int i = 0; i < tags.size(); ++i) {

      }

      // TODO process all sequences after last tag

      // We produced all input sequences except the last one, because it has no successor.
      return noutput_items-2;
    }

  } /* namespace digital */
} /* namespace gr */


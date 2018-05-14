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

#include <gnuradio/io_signature.h>
#include "alamouti_encoder_cc_impl.h"

using namespace boost;

namespace gr {
  namespace digital {

    alamouti_encoder_cc::sptr
    alamouti_encoder_cc::make()
    {
      return gnuradio::get_initial_sptr
        (new alamouti_encoder_cc_impl());
    }

    /*
     * The private constructor
     */
    alamouti_encoder_cc_impl::alamouti_encoder_cc_impl()
      : gr::sync_block("alamouti_encoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(2, 2, sizeof(gr_complex)))
    {
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2);
    }

    /*
     * Our virtual destructor.
     */
    alamouti_encoder_cc_impl::~alamouti_encoder_cc_impl()
    {
    }

    int
    alamouti_encoder_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out1 = (gr_complex *) output_items[0];
      gr_complex *out2 = (gr_complex *) output_items[1];

      // Copy data to first period of each transmission sequence.
      for (int i = 0; i < noutput_items; i+=2) {
        out1[i] = in[i];
        out2[i] = in[i+1];
      }
      // Write conjugated (and for branch 2 negated) data to 2. period of each transmission sequence.
      for (int i = 1; i < noutput_items; i+=2) {
        out1[i] = -std::conj(in[i]);
        out2[i] = std::conj(in[i-1]);
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


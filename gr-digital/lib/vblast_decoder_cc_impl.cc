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
#include "vblast_decoder_cc_impl.h"

namespace gr {
  namespace digital {

    const std::string vblast_decoder_cc_impl::s = "csi";
    const pmt::pmt_t vblast_decoder_cc_impl::d_key = pmt::string_to_symbol(s);

    vblast_decoder_cc::sptr
    vblast_decoder_cc::make(uint16_t num_inputs)
    {
      return gnuradio::get_initial_sptr
        (new vblast_decoder_cc_impl(num_inputs));
    }

    /*
     * The private constructor
     */
    vblast_decoder_cc_impl::vblast_decoder_cc_impl(uint16_t num_inputs)
      : gr::sync_interpolator("vblast_decoder_cc",
              gr::io_signature::make(num_inputs, num_inputs, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), num_inputs),
        d_num_inputs(num_inputs)
    {}

    /*
     * Our virtual destructor.
     */
    vblast_decoder_cc_impl::~vblast_decoder_cc_impl()
    {
    }

    int
    vblast_decoder_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];
      uint16_t nprocessed = 0; // Number of read and written items.

      // Collect all tags of the input buffer with key "csi" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);

      uint16_t symbol_length; // Number of items in the current symbol.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items;
        //TODO Process symbol.
        nprocessed += symbol_length;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
           * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          //TODO Process_symbol.
          nprocessed += symbol_length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            symbol_length = tags[i + 1].offset - tags[i].offset;
          } else {
            symbol_length = noutput_items - tags[i].offset + nitems_read(0);
          }
          // Get CSI from tag.
          d_csi = pmt::c32vector_elements(tags[i].value);
          // Calculate the weighting vector for the next symbol with the received CSI.
          //TODO Calculate the weighting vector.
          // Process the symbol with the calculated weighting vector.
          //TODO Process_symbol.
          nprocessed += symbol_length;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


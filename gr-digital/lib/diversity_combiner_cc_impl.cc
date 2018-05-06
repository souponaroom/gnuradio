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

#include <boost/format.hpp>
#include <gnuradio/io_signature.h>
#include "diversity_combiner_cc_impl.h"
#include <volk/volk.h>
#include <pmt/pmt.h>

using namespace boost;

namespace gr {
  namespace digital {

    diversity_combiner_cc::sptr
    diversity_combiner_cc::make(uint16_t num_inputs, uint16_t vlen, uint8_t combining_technique)
    {
      return gnuradio::get_initial_sptr
        (new diversity_combiner_cc_impl(num_inputs, vlen, combining_technique));
    }

    /*
     * The private constructor
     */
    diversity_combiner_cc_impl::diversity_combiner_cc_impl(uint16_t num_inputs, uint16_t vlen, uint8_t combining_technique)
      : gr::sync_block("diversity_combiner_cc",
              gr::io_signature::make(num_inputs, num_inputs, vlen*sizeof(gr_complex)),
              gr::io_signature::make(1, 1, vlen*sizeof(gr_complex))),
        d_num_inputs(num_inputs),
        d_vlen(vlen),
        d_combining_technique(combining_technique),
        d_best_path(0)
    {
      // Initialize csi vector.
      for (int i = 0; i < num_inputs; ++i) {
        d_csi.push_back(1.0);
        d_csi_squared.push_back(1.0);
      }
      // Set tag propagation policy to 'All to All'.
      set_tag_propagation_policy(TPP_ALL_TO_ALL);
      // Allocate space for repeating volk operations.
      //unsigned int alignment = volk_get_alignment();
      //d_csi_squared = (float *) volk_malloc(sizeof(float) * num_inputs, alignment);
    }

    /*
     * Our virtual destructor.
     */
    diversity_combiner_cc_impl::~diversity_combiner_cc_impl()
    {
    }

    void
    diversity_combiner_cc_impl::process_symbol(gr_vector_const_void_star input, gr_complex* out, uint16_t offset, uint16_t length){
      switch (d_combining_technique) {
        case 0: { // Selection combining
          // Search for path coefficient with maximal magnitude.
          // Therefore we calculate the magnitude square out of the complex CSI vector.
          for (int i = 0; i < d_num_inputs; ++i) {
            d_csi_squared[i] = std::abs(d_csi[i]);
          }
          // Select maximum value.
          d_best_path = std::distance(d_csi_squared.begin(),
                                             std::max_element(d_csi_squared.begin(),
                                                              d_csi_squared.end()));
          const gr_complex *in = &((const gr_complex *) input[d_best_path])[offset];
          // Copy items of the current symbol from best_path to output.
          memcpy(out, in, length * sizeof(gr_complex) * d_vlen);
          break;
        }
        case 1: { // Maximum-Ration combining
          break;
        }
      }
    }

    int
    diversity_combiner_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      // Define output pointer 'out' and its control variable 'nprocessed'.
      gr_complex *out = (gr_complex *) output_items[0];
      uint16_t nprocessed = 0;

      // Collect all tags of the input buffer with key "csi" in the vector 'tags'.
      std::vector <gr::tag_t> tags;
      const std::string s = "csi";
      pmt::pmt_t d_key = pmt::string_to_symbol(s);
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);

      uint16_t symbol_length;

      // Handle samples before the first tag (with memory).
      if(tags.size() == 0){ // Input buffer includes no tags at all.
        symbol_length = noutput_items;
        switch (d_combining_technique) {
          case 0: { // Selection combining
            const gr_complex *in = &((const gr_complex *) input_items[d_best_path])[nprocessed];
            // Copy items of the current symbol from best_path to output.
            memcpy(&out[nprocessed], in, symbol_length * sizeof(gr_complex) * d_vlen);
            nprocessed += symbol_length;
            break;
          }
          case 1: { // Maximum-Ratio Combining
            break;
          }
        }
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
           * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          switch (d_combining_technique) {
            case 0: { // Selection combining
              const gr_complex *in = &((const gr_complex *) input_items[d_best_path])[nprocessed];
              // Copy items of the current symbol from best_path to output.
              memcpy(&out[nprocessed], in, symbol_length * sizeof(gr_complex) * d_vlen);
              nprocessed += symbol_length;
              break;
            }
            case 1: {
              break;
            }
          }
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
          // Process the next symbol with the received CSI.
          process_symbol(input_items, &out[nprocessed], nprocessed, symbol_length);
          nprocessed += symbol_length;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


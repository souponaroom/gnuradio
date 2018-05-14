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
#include <stdio.h>
#include <string>
#include <math.h>
#include <numeric>
#include <complex>

using namespace boost;

namespace gr {
  namespace digital {

    const std::string diversity_combiner_cc_impl::s = "csi";
    const pmt::pmt_t diversity_combiner_cc_impl::d_key = pmt::string_to_symbol(s);

    diversity_combiner_cc::sptr
    diversity_combiner_cc::make(uint16_t num_inputs, uint16_t vlen, std::string combining_technique)
    {
      return gnuradio::get_initial_sptr
        (new diversity_combiner_cc_impl(num_inputs, vlen, combining_technique));
    }

    /*
     * The private constructor
     */
    diversity_combiner_cc_impl::diversity_combiner_cc_impl(uint16_t num_inputs,
                                                           uint16_t vlen,
                                                           std::string combining_technique)
      : gr::sync_block("diversity_combiner_cc",
              gr::io_signature::make(num_inputs, num_inputs, vlen*sizeof(gr_complex)),
              gr::io_signature::make(1, 1, vlen*sizeof(gr_complex))),
        d_num_inputs(num_inputs),
        d_vlen(vlen),
        d_combining_technique(combining_technique),
        d_best_path(0) // Initially, the SC algorithm selects channel 0, until it is changed by CSI.
    {
      // Resize vectors csi, csi_squared and mrc_weighting.
      d_csi.resize(num_inputs);
      d_csi_squared.resize(num_inputs);
      d_mrc_weighting.resize(num_inputs);
      // Initially, the MRC algorithm weights each channel equally, until it is changed by CSI.
      std::fill(d_mrc_weighting.begin(), d_mrc_weighting.end(), 1.0/num_inputs);
      // Set tag propagation policy to 'All to All'.
      set_tag_propagation_policy(TPP_ALL_TO_ALL);
    }

    /*
     * Our virtual destructor.
     */
    diversity_combiner_cc_impl::~diversity_combiner_cc_impl()
    {
    }

    void
    diversity_combiner_cc_impl::combine_inputs(gr_vector_const_void_star input,
                                               gr_complex* out,
                                               uint64_t offset,
                                               uint64_t length){
      // Calculate the magnitude square of the complex CSI vector.
      for (int i = 0; i < d_num_inputs; ++i) {
        d_csi_squared[i] = std::norm(d_csi[i]);
      }
      if(d_combining_technique.compare("SC") == 0) { // Selection combining
        // Search for path coefficient with maximum magnitude.
        d_best_path = std::distance(d_csi_squared.begin(),
                                    std::max_element(d_csi_squared.begin(),
                                                     d_csi_squared.end()));
      }
      else if(d_combining_technique.compare("MRC") == 0) { // Maximum-Ratio combining
        // Calculate the normalized weighting coefficients.
        float total_path_energy = std::accumulate(d_csi_squared.begin(), d_csi_squared.end(), 0.0);
        for (int i = 0; i < d_num_inputs; ++i) {
          d_mrc_weighting[i] = std::polar(std::sqrt(d_csi_squared[i]/total_path_energy), -std::arg(d_csi[i]));
        }
      }

    }

    void
    diversity_combiner_cc_impl::process_symbol(gr_vector_const_void_star input,
                                               gr_complex* out,
                                               uint64_t offset,
                                               uint64_t length){
      if(d_combining_technique.compare("SC") == 0) { // Selection combining
        const gr_complex *in = &((const gr_complex *) input[d_best_path])[offset];
        // Copy items of the current symbol from best_path to output.
        std::copy(in, &in[length*d_vlen], out);
      }
      else if(d_combining_technique.compare("MRC") == 0) { // Maximum-Ratio combining
        // Calculate the output stream as the weighted sum of the input streams.
        std::fill(out, &out[length*d_vlen], 0.0);
        for (int inport = 0; inport < d_num_inputs; ++inport) {
          const gr_complex *in = &((const gr_complex *) input[inport])[offset];
          for (unsigned int l = 0; l < length * d_vlen; ++l) {
            out[l] += d_mrc_weighting[inport] * in[l];
          }
        }
      }
    }

    int
    diversity_combiner_cc_impl::work(int noutput_items,
                                     gr_vector_const_void_star &input_items,
                                     gr_vector_void_star &output_items)
    {
      gr_complex *out = (gr_complex *) output_items[0];
      uint16_t nprocessed = 0; // Number of read and written items (vectors, if d_vlen > 1).

      // Collect all tags of the input buffer with key "csi" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);

      uint16_t symbol_length; // Number of items in the current symbol.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items;
        process_symbol(input_items, out, 0, symbol_length);
        nprocessed += symbol_length;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
           * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          process_symbol(input_items, out, 0, symbol_length);
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
          combine_inputs(input_items, &out[nprocessed*d_vlen], nprocessed*d_vlen, symbol_length);
          // Process the symbol with the calculated weighting vector.
          process_symbol(input_items, &out[nprocessed*d_vlen], nprocessed*d_vlen, symbol_length);
          nprocessed += symbol_length;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


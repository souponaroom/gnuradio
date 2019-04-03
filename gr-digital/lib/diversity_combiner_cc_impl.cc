/* -*- c++ -*- */
/* 
 * Copyright 2018, 2019 Free Software Foundation, Inc.
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

    diversity_combiner_cc::sptr
    diversity_combiner_cc::make(uint16_t num_inputs, uint16_t vlen,
                                std::string combining_technique,
                                const std::string &csi_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new diversity_combiner_cc_impl(num_inputs, vlen, combining_technique, csi_tag_key));
    }

    /*
     * The private constructor
     */
    diversity_combiner_cc_impl::diversity_combiner_cc_impl(uint16_t num_inputs,
                                                           uint16_t vlen,
                                                           std::string combining_technique,
                                                           const std::string &csi_tag_key)
      : gr::sync_interpolator("diversity_combiner_cc",
              gr::io_signature::make(num_inputs, num_inputs, vlen*sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), vlen),
        d_num_inputs(num_inputs),
        d_vlen(vlen),
        d_combining_technique(combining_technique),
        d_csi_key(pmt::string_to_symbol(csi_tag_key))
    {
      d_csi = std::vector<std::vector<std::vector<gr_complex> > >(vlen, std::vector<std::vector<gr_complex> >(num_inputs, std::vector<gr_complex> (1, 1.0)));
      d_csi_squared = std::vector<std::vector<float> >(vlen, std::vector<float>(num_inputs, 1.0));
      d_best_path = std::vector<uint16_t> (vlen, 0); // Initially, the SC algorithm selects channel 0, until it is changed by CSI.
      // Initially, the MRC algorithm weights each channel equally, until it is changed by CSI.
      d_mrc_weighting = std::vector<std::vector<gr_complex> >(vlen, std::vector<gr_complex>(num_inputs, 1.0/num_inputs));
      // Set tag propagation policy to 'All to All'.
      set_tag_propagation_policy(TPP_DONT);
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
      for (int k = 0; k < d_vlen; ++k) {
        for (int i = 0; i < d_num_inputs; ++i) {
          d_csi_squared[k][i] = std::norm(d_csi[k][i][0]);
        }
      }
      if(d_combining_technique.compare("SC") == 0) { // Selection combining
        // Search for path coefficient with maximum magnitude.
        for (int k = 0; k < d_vlen; ++k) {
          d_best_path[k] = std::distance(d_csi_squared[k].begin(),
                                      std::max_element(d_csi_squared[k].begin(),
                                                       d_csi_squared[k].end()));
        }
      }
      else if(d_combining_technique.compare("MRC") == 0) { // Maximum-Ratio combining
        // Calculate the normalized weighting coefficients.
        for (int k = 0; k < d_vlen; ++k){
          float total_path_energy = std::accumulate(d_csi_squared[k].begin(), d_csi_squared[k].end(), 0.0);
          for (int i = 0; i < d_num_inputs; ++i) {
            d_mrc_weighting[k][i] = std::polar(std::sqrt(d_csi_squared[k][i]/total_path_energy), -std::arg(d_csi[k][i][0]));
          }
        }
      }

    }

    void
    diversity_combiner_cc_impl::process_symbol(gr_vector_const_void_star input,
                                               gr_complex* out,
                                               uint64_t offset,
                                               uint64_t length){
      if(d_combining_technique.compare("SC") == 0) { // Selection combining
        // Copy items of the current symbol from best_path to output.
        for (int k = 0; k < d_vlen; ++k) {
          const gr_complex *in = &((const gr_complex *) input[d_best_path[k]])[offset];
          for (int i = 0; i < length; ++i) {
            out[d_vlen*i+k] = in[d_vlen*i+k]/d_csi[k][d_best_path[k]][0];
          }
        }
      }
      else if(d_combining_technique.compare("MRC") == 0) { // Maximum-Ratio combining
        // Calculate the output stream as the weighted sum of the input streams.
        std::fill(out, &out[length * d_vlen], 0.0);
        for (int k = 0; k < d_vlen; ++k) {
          for (int inport = 0; inport < d_num_inputs; ++inport) {
            const gr_complex *in = &((const gr_complex *) input[inport])[offset];
            for (unsigned int l = 0; l < length; ++l) {
              out[l*d_vlen+k] += d_mrc_weighting[k][inport] * in[l*d_vlen+k];
            }
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
      get_tags_in_window(d_tags, 0, 0, noutput_items/d_vlen, d_csi_key);

      uint16_t symbol_length; // Number of items in the current symbol.

      if(d_tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items/d_vlen;
        process_symbol(input_items, out, 0, symbol_length);
        nprocessed += symbol_length;
      } else { // Input buffer includes tags.
        if (d_tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
           * which belong to the previous symbol. */
          symbol_length = d_tags[0].offset - nitems_read(0);
          process_symbol(input_items, out, 0, symbol_length);
          nprocessed += symbol_length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < d_tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < d_tags.size() - 1) {
            symbol_length = d_tags[i + 1].offset - d_tags[i].offset;
          } else {
            symbol_length = noutput_items/d_vlen - d_tags[i].offset + nitems_read(0);
          }
          // Get CSI from tag.
          for (unsigned int k = 0; k < pmt::length(d_tags[i].value); ++k) {
            pmt::pmt_t carrier_csi = pmt::vector_ref(d_tags[i].value, k);
            for (unsigned int j = 0; j < pmt::length(carrier_csi); ++j) {
              d_csi[k][j] = pmt::c32vector_elements(pmt::vector_ref(carrier_csi, j));
            }
          }
          // Calculate the weighting vector for the next symbol with the received CSI.
          combine_inputs(input_items, &out[nprocessed*d_vlen], nprocessed*d_vlen, symbol_length);
          // Process the symbol with the calculated weighting vector.
          process_symbol(input_items, &out[nprocessed*d_vlen], nprocessed*d_vlen, symbol_length);
          nprocessed += symbol_length;
        }
      }
      // Propagate all other tags (except the CSI tags which were changed) manually.
      std::vector <gr::tag_t> tags;
      get_tags_in_window(tags, 0, 0, noutput_items/d_vlen);
      for (int l = 0; l < tags.size(); ++l) {
        if (tags[l].key != d_csi_key) {
          add_item_tag(0, tags[l].offset*d_vlen, tags[l].key, tags[l].value);
        }
      }
      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


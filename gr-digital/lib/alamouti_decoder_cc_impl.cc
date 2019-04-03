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
#include "alamouti_decoder_cc_impl.h"

using namespace boost;

namespace gr {
  namespace digital {

    alamouti_decoder_cc::sptr
    alamouti_decoder_cc::make(uint32_t vlen, const std::string &csi_tag_key)
    {
      return gnuradio::get_initial_sptr
        (new alamouti_decoder_cc_impl(vlen, csi_tag_key));
    }

    /*
     * The private constructor
     */
    alamouti_decoder_cc_impl::alamouti_decoder_cc_impl(uint32_t vlen,
                                                       const std::string &csi_tag_key)
      : gr::sync_interpolator("alamouti_decoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)*vlen),
              gr::io_signature::make(1, 1, sizeof(gr_complex)), vlen),
        d_vlen(vlen),
        d_csi_key(pmt::string_to_symbol(csi_tag_key))
    {
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2*vlen);
      // Check if vlen is an even number.
      if (vlen != 1 && vlen%2 != 0){
        throw std::invalid_argument("Vector length must be an even number.");
      }
      // Init CSI array.
      d_csi = std::vector<std::vector<std::vector<gr_complex> > >(vlen, std::vector<std::vector<gr_complex> >(1, std::vector<gr_complex> (2, 1.0)));
      // Set tag propagation policy to 'All to All'.
      set_tag_propagation_policy(TPP_DONT);
    }

    /*
     * Our virtual destructor.
     */
    alamouti_decoder_cc_impl::~alamouti_decoder_cc_impl()
    {
    }

    void
    alamouti_decoder_cc_impl::decode_symbol(const gr_complex* in,
                                            gr_complex* out,
                                            uint32_t length){
      // Iterate over received sequences (= 2 complex symbols).
      for (unsigned int i = 0; i < length; i+=2) {
        for (int k = 0; k < d_vlen; ++k) {
          // Calculate the sum of the energy of both branches.
          float total_branch_energy = std::norm(d_csi[k][0][0]) + std::norm(d_csi[k][0][1]);
          // Calculate an estimation for the transmission sequence.
          out[i*d_vlen+k] = (std::conj(d_csi[k][0][0])*in[i*d_vlen+k]+d_csi[k][0][1]*std::conj(in[(i+1)*d_vlen+k]))/total_branch_energy;
          out[(i+1)*d_vlen+k] = (std::conj(d_csi[k][0][1])*in[i*d_vlen+k]-d_csi[k][0][0]*std::conj(in[(i+1)*d_vlen+k]))/total_branch_energy;
        }
      }
    }

    int
    alamouti_decoder_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      uint32_t nprocessed = 0; // Number of read and written items.

      // Collect all tags of the input buffer with key "csi" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items/d_vlen, d_csi_key);

      uint32_t symbol_length; // Number of items in the current symbol.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items/d_vlen;
        decode_symbol(in, out, symbol_length);
        nprocessed += symbol_length*d_vlen;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
+          * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          // Check if the next tag is on an uneven position.
          if(tags[0].offset%2 != 0){
            // The CSI is updated with the start of the next sequence (=next even sample).
            ++symbol_length;
          }
          decode_symbol(in, out, symbol_length);
          nprocessed += symbol_length*d_vlen;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag in the buffer.
            // Check if there are 2 tags on the same item.
            if(tags[i+1].offset <= tags[i].offset){
              continue;
            }
            symbol_length = tags[i + 1].offset - nitems_read(0) - nprocessed/d_vlen;
            // Check if the next tag is on an uneven position.
            if(symbol_length%2 != 0){
              // The CSI is updated with the start of the next sequence (=next even sample).
              ++symbol_length;
            }
          } else {
            // This is the last tag in the buffer.
            symbol_length = (noutput_items - nprocessed)/d_vlen;
          }
          // Get CSI from 'csi' tag.
          for (unsigned int k = 0; k < pmt::length(tags[i].value); ++k) {
            pmt::pmt_t carrier_csi = pmt::vector_ref(tags[i].value, k);
            for (unsigned int j = 0; j < pmt::length(carrier_csi); ++j) {
              d_csi[k][j] = pmt::c32vector_elements(pmt::vector_ref(carrier_csi, j));
            }
          }
          // Process the symbol with the calculated weighting vector.
          decode_symbol(&in[nprocessed], &out[nprocessed], symbol_length);
          nprocessed += symbol_length*d_vlen;
        }
      }
      // Propagate all other tags (except the CSI tags which were changed) manually.
      std::vector <gr::tag_t> tags;
      get_tags_in_window(tags, 0, 0, noutput_items);
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


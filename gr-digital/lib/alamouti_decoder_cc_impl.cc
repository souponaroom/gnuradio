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
#include "alamouti_decoder_cc_impl.h"

using namespace boost;

namespace gr {
  namespace digital {

    const std::string alamouti_decoder_cc_impl::s = "csi";
    const pmt::pmt_t alamouti_decoder_cc_impl::d_key = pmt::string_to_symbol(s);

    alamouti_decoder_cc::sptr
    alamouti_decoder_cc::make()
    {
      return gnuradio::get_initial_sptr
        (new alamouti_decoder_cc_impl());
    }

    /*
     * The private constructor
     */
    alamouti_decoder_cc_impl::alamouti_decoder_cc_impl()
      : gr::sync_block("alamouti_decoder_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
      /* Set the number of input and output items to a multiple of 2,
       * because the Alamouti algorithm processes sequences of 2 complex symbols.
       */
      set_output_multiple(2);
      // Init CSI array.
      d_csi.push_back(1.0);
      d_csi.push_back(1.0);
      // Set tag propagation policy to 'All to All'.
      set_tag_propagation_policy(TPP_ALL_TO_ALL);
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
        // Calculate the sum of the energy of both branches.
        float total_branch_energy = std::norm(d_csi[0]) + std::norm(d_csi[1]);
        // Calculate an estimation for the transmission sequence.
        out[i] = (std::conj(d_csi[0])*in[i] + d_csi[1]*std::conj(in[i+1]))/total_branch_energy;
        out[i+1] = (std::conj(d_csi[1])*in[i] - d_csi[0]*std::conj(in[i+1]))/total_branch_energy;
      }
    }

    int
    alamouti_decoder_cc_impl::work(int noutput_items,
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
        decode_symbol(in, out, symbol_length);
        nprocessed += symbol_length;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
+          * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          // Check if the next tag is on an uneven position.
          if(tags[0].offset%2 != 0){
            // This should be prevented by the system developer in most cases.
            GR_LOG_WARN(d_logger, format("Detected \'csi\' tag on uneven position (tag[%d].offset = %d).\n "
                                         "The Alamouti scheme works on sequences of 2 samples. "
                                         "If you are not really sure what you are doing, "
                                         "you should only set 'csi' tags on even sample positions.")
                                  %0 %tags[0].offset);
            // The CSI is updated with the start of the next sequence (=next even sample).
            ++symbol_length;
          }
          decode_symbol(in, out, symbol_length);
          nprocessed += symbol_length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag in the buffer.
            symbol_length = tags[i + 1].offset - nitems_read(0) - nprocessed;
            // Check if the next tag is on an uneven position (which it should usually not).
            if(symbol_length%2 != 0){
              // This should be prevented by the system developer in most cases.
              GR_LOG_WARN(d_logger, format("Detected \'csi\' tag on uneven position (tag[%d].offset = %d). \n"
                                                   "The Alamouti scheme works on sequences of 2 samples. "
                                                   "If you are not really sure what you are doing, "
                                                   "you should only set 'csi' tags on even sample positions.")
                                    %i %tags[i].offset);
              // The CSI is updated with the start of the next sequence (=next even sample).
              ++symbol_length;
            }
          } else {
            // This is the last tag in the buffer.
            symbol_length = noutput_items - nprocessed;
          }
          // Get CSI from tag.
          d_csi = pmt::c32vector_elements(tags[i].value);
          // Process the symbol with the calculated weighting vector.
          decode_symbol(&in[nprocessed], &out[nprocessed], symbol_length);
          nprocessed += symbol_length;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


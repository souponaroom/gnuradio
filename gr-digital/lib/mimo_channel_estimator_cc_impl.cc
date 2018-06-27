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
#include "mimo_channel_estimator_cc_impl.h"

#define _USE_MATH_DEFINES
#include <math.h>

#if (EIGEN3_ENABLED)
#include <eigen3/Eigen/Dense>
#endif

namespace gr {
  namespace digital {

    const pmt::pmt_t mimo_channel_estimator_cc_impl::d_key = pmt::string_to_symbol("pilot");

    mimo_channel_estimator_cc::sptr
    mimo_channel_estimator_cc::make(uint16_t num_inputs, std::vector<std::vector<gr_complex> > training_sequence)
    {
      return gnuradio::get_initial_sptr
        (new mimo_channel_estimator_cc_impl(num_inputs, training_sequence));
    }

    /*
     * The private constructor
     */
    mimo_channel_estimator_cc_impl::mimo_channel_estimator_cc_impl(uint16_t num_inputs, std::vector<std::vector<gr_complex> > training_sequence)
      : gr::block("mimo_channel_estimator_cc",
              gr::io_signature::make(num_inputs, num_inputs, sizeof(gr_complex)),
              gr::io_signature::make(num_inputs, num_inputs, sizeof(gr_complex))),
        d_num_inputs(num_inputs),
        d_training_sequence(training_sequence)
    {
      d_training_length = d_training_sequence[0].size();
    }

    /*
     * Our virtual destructor.
     */
    mimo_channel_estimator_cc_impl::~mimo_channel_estimator_cc_impl()
    {
    }

    void
    mimo_channel_estimator_cc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    void
    mimo_channel_estimator_cc_impl::copy_symbols(gr_vector_const_void_star &input_items,
                                                 gr_vector_void_star &output_items,
                                                 uint32_t symbol_length,
                                                 uint32_t reading_offset,
                                                 uint32_t writing_offset) {
      for (int i = 0; i < d_num_inputs; ++i) {
        const gr_complex *in = (const gr_complex *) input_items[i];
        gr_complex *out = (gr_complex *) output_items[i];
        memcpy(&out[writing_offset], &in[reading_offset], symbol_length);
      }
    }

    void
    mimo_channel_estimator_cc_impl::estimate_channel(gr_vector_const_void_star &input_items,
                                                     uint32_t reading_offset) {
      switch (d_num_inputs) {
        case 1: {
          // Init CSI vector.
          std::vector<std::vector<gr_complex> > csi (1, std::vector<gr_complex> (1, 0.0));
          // Fill CSI vector with precalculated ML channel estimation.
          for (int i = 0; i < d_training_length; ++i) {
            csi[0][0] += ((const gr_complex *) input_items[0])[i];
          }
          csi[0][0] *= 1./ d_training_length;
          break;
        }
        case 2: {
          // Init CSI vector.
          std::vector<std::vector<gr_complex> > csi (2, std::vector<gr_complex> (2, 0.0));
          // Fill CSI vector with precalculated ML channel estimation.
          for (int i = 0; i < d_training_length; ++i) {
            csi[0][0] += ((const gr_complex *) input_items[0])[i];
            csi[0][1] += ((const gr_complex *) input_items[0])[i] * (gr_complex) std::polar(1.0, 2*M_PI*i / d_training_length);
            csi[1][0] += ((const gr_complex *) input_items[1])[i];
            csi[1][1] += ((const gr_complex *) input_items[1])[i] * (gr_complex) std::polar(1.0, 2*M_PI*i / d_training_length);
          }
          // Multiply with factor.
          // TODO Integrate SNR in factor!!!
          csi[0][0] *= M_SQRT2 / d_training_length;
          csi[0][1] *= M_SQRT2 / d_training_length;
          csi[1][0] *= M_SQRT2 / d_training_length;
          csi[1][1] *= M_SQRT2 / d_training_length;
          break;
        }
        default: {

        }
      }
    }

    int
    mimo_channel_estimator_cc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      // Collect all tags of the input buffer with key "csi" in the vector 'tags'.
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);

      uint32_t nconsumed = 0; // Number of read items.
      uint32_t nwritten = 0; // Number of written items.
      uint32_t symbol_length; // Number of items in the current symbol.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Handle all samples in buffer as they belong to the current symbol.
        symbol_length = noutput_items;
        copy_symbols(input_items, output_items, symbol_length, 0, 0);
        nconsumed = noutput_items;
        nwritten = noutput_items;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0){
          /* There are items in the input buffer, before the first tag arrives,
+          * which belong to the previous symbol. */
          symbol_length = tags[0].offset - nitems_read(0);
          copy_symbols(input_items, output_items, symbol_length, 0, 0);
          nconsumed += symbol_length;
          nwritten += symbol_length - d_training_length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i) {
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag in the buffer.
            symbol_length = tags[i + 1].offset - nitems_read(0) - nwritten;
          } else {
            // This is the last tag in the buffer.
            symbol_length = noutput_items - nwritten;
          }
          // Copy symbols to output.
          copy_symbols(input_items, output_items, symbol_length, nconsumed, nwritten);
          // Estimate MIMO channel and write CSI tag to stream.
          // TODO handel other previous training symbols (Schmidl & Cox) via pilot offset or dumping
          estimate_channel(input_items, nconsumed);

          nconsumed += symbol_length;
          nwritten += symbol_length - d_training_length;
        }
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


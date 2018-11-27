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

#define _USE_MATH_DEFINES
#define M_TWOPI (2*M_PI)
#include <cmath>
#include <gnuradio/expj.h>
#include <gnuradio/io_signature.h>
#include "ofdm_correct_carrier_freq_offset_vcvc_impl.h"

namespace gr {
  namespace digital {

    ofdm_correct_carrier_freq_offset_vcvc::sptr
    ofdm_correct_carrier_freq_offset_vcvc::make(uint16_t n, uint32_t fft_len, uint32_t cp_len, std::string carrier_freq_offset_key)
    {
      return gnuradio::get_initial_sptr
        (new ofdm_correct_carrier_freq_offset_vcvc_impl(n, fft_len, cp_len, carrier_freq_offset_key));
    }

    /*
     * The private constructor
     */
    ofdm_correct_carrier_freq_offset_vcvc_impl::ofdm_correct_carrier_freq_offset_vcvc_impl(
            uint16_t n,
            uint32_t fft_len,
            uint32_t cp_len,
            std::string carrier_freq_offset_key)
      : gr::sync_block("ofdm_correct_carrier_freq_offset_vcvc",
              gr::io_signature::make(n, n, sizeof(gr_complex)*fft_len),
              gr::io_signature::make(n, n, sizeof(gr_complex)*fft_len)),
        d_n(n),
        d_fft_len(fft_len),
        d_cp_len(cp_len),
        d_carrier_freq_offset_key(carrier_freq_offset_key),
        d_key(pmt::string_to_symbol(carrier_freq_offset_key)),
        d_carrier_offset(0)
    {
      // Set tag propagation policy.
      set_tag_propagation_policy(TPP_ONE_TO_ONE);
    }

    /*
     * Our virtual destructor.
     */
    ofdm_correct_carrier_freq_offset_vcvc_impl::~ofdm_correct_carrier_freq_offset_vcvc_impl()
    {
    }

    void
    ofdm_correct_carrier_freq_offset_vcvc_impl::correct_offset(gr_vector_const_void_star &input_items,
                                                               gr_vector_void_star &output_items,
                                                               uint32_t offset, uint32_t length) {
      // Iterate over branches.
      for (int n = 0; n < d_n; ++n) {
        const gr_complex *in = &(((const gr_complex *) input_items[n])[offset]);
        gr_complex *out = &(((gr_complex *) output_items[n])[offset]);

        // Copy the fft vectors such that the symbols are shifted to the correct position.
        if (d_carrier_offset < 0) {
          memset((void *) out, 0x00, sizeof(gr_complex) * (-d_carrier_offset));
          memcpy((void *) &out[-d_carrier_offset], (void *) in,
                 sizeof(gr_complex) * (d_fft_len * length + d_carrier_offset));
        } else {
          memset((void *) (out + d_fft_len * length - d_carrier_offset),
                 0x00, sizeof(gr_complex) * d_carrier_offset);
          memcpy((void *) out, (void *) (in + d_carrier_offset),
                 sizeof(gr_complex) * (d_fft_len * length - d_carrier_offset));
        }
        /* The cyclic prefix was cut out somewhere before this block.
         * But the carrier frequency offset was not yet corrected back then.
         * (Because it is corrected right now). This leads to a phase shift
         * between the FFT vectors. */
        // Correct this phase shift.
        gr_complex phase_correction;
        for (unsigned int i = 0; i < length; i++) {
          phase_correction = gr_expj(-M_TWOPI * d_carrier_offset * d_cp_len / d_fft_len * (i + 1));
          for (unsigned int k = 0; k < d_fft_len; k++) {
            out[i * d_fft_len + k] *= phase_correction;
          }
        }
      }
    }

    int
    ofdm_correct_carrier_freq_offset_vcvc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      // Collect all tags of the input buffer in the vector 'tags'.
      std::vector <gr::tag_t> tags;
      get_tags_in_window(tags, 0, 0, noutput_items, d_key);
      uint32_t length; // Number of items until the next tag arrives.
      uint32_t nprocessed = 0; // Number of read and written items.

      if(tags.size() == 0){ // Input buffer includes no tags at all.
        // Correct int freq offset of the whole buffer.
        correct_offset(input_items, output_items, 0, noutput_items);
        nprocessed = noutput_items;
      } else { // Input buffer includes tags.
        if (tags[0].offset - nitems_read(0) > 0) {
          /* There are items in the input buffer, before the first tag arrives,
           * which belong to the previous frame. */
          length = (tags[0].offset - nitems_read(0));
          // Correct int freq offset.
          correct_offset(input_items, output_items, 0, length);
          nprocessed = length;
        }
        // Iterate over tags in buffer.
        for (unsigned int i = 0; i < tags.size(); ++i){
          // Read the new carrier offset.
          d_carrier_offset = pmt::to_long(tags[i].value);
          // Calculate the number of items before the next tag.
          if (i < tags.size() - 1) {
            // This is not the last tag.
            length = (tags[i + 1].offset - tags[i].offset);
          } else {
            // This is the last tag.
            length = noutput_items - (tags[i].offset - nitems_read(0));
          }
          // Correct int freq offset.
          correct_offset(input_items, output_items, nprocessed*d_fft_len, length);
          nprocessed += length;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


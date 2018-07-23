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
#include <cmath>

#include <gnuradio/io_signature.h>
#include "mimo_ofdm_synchronizer_fbcvc_impl.h"

namespace gr {
  namespace digital {

    mimo_ofdm_synchronizer_fbcvc::sptr
    mimo_ofdm_synchronizer_fbcvc::make(const std::vector<gr_complex> &sync_symbol1,
                                       const std::vector<gr_complex> &sync_symbol2,
                                       uint16_t cp_len,
                                       uint16_t n)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_synchronizer_fbcvc_impl(sync_symbol1,
                                               sync_symbol2,
                                               cp_len,
                                               n));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_synchronizer_fbcvc_impl::mimo_ofdm_synchronizer_fbcvc_impl(
            const std::vector<gr_complex> &sync_symbol1,
            const std::vector<gr_complex> &sync_symbol2,
            uint16_t cp_len,
            uint16_t n)
      : gr::block("mimo_ofdm_synchronizer_fbcvc",
              gr::io_signature::make3(3+n, 3+n, sizeof(float), sizeof(unsigned char), sizeof(gr_complex)),
              gr::io_signature::make(n, n, sizeof(gr_complex) * sync_symbol1.size())),
        d_cp_len(cp_len),
        d_n(n),
        d_fft_len(sync_symbol1.size()),
        d_sync_sym_count(0),
        d_phase(0.)
    {
      d_symbol_len = d_fft_len + d_cp_len;
    }

    /*
     * Our virtual destructor.
     */
    mimo_ofdm_synchronizer_fbcvc_impl::~mimo_ofdm_synchronizer_fbcvc_impl()
    {
    }

    void
    mimo_ofdm_synchronizer_fbcvc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      for (int i = 0; i < 3+d_n; ++i) {
        ninput_items_required[i] = noutput_items*(d_fft_len+d_cp_len);
      }
    }

    void
    mimo_ofdm_synchronizer_fbcvc_impl::rotate_phase(const float *fine_freq_off,
                                                    uint16_t rotation_length) {
      for (int i = 0; i < rotation_length; ++i) {
        d_phase += fine_freq_off[i] * 2.0 / d_fft_len;
      }
      // Place the phase in [0, pi)
      d_phase = std::fmod(d_phase, 2.*M_PI);
    }

    int
    mimo_ofdm_synchronizer_fbcvc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const float *fine_freq_off = (const float *) input_items[0];
      const unsigned char *trigger = (const unsigned char *) input_items[1];
      const gr_complex *ref_sig = (const gr_complex *) input_items[2];

      uint32_t nconsumed = 0;
      uint32_t nwritten = 0;

      while (nconsumed < noutput_items*(d_symbol_len)) {
        // Check if the next symbol start hast a trigger.

        // Search for a trigger in following items of the symbols (except the first one).
        uint16_t next_trigger = 0;
        for (unsigned int k = 1; k < d_fft_len; ++k) {
          if (trigger[nconsumed + k] == 1) {
            next_trigger = k;
            // We found a trigger in the symbol. Apparently, we are out of sync.
            if (trigger[nconsumed] == 1){
              /* There is also a trigger at the beginning of the symbol. We would
               * dumpy the first 2 symbols (sync symbols), but due the trigger
               * later in the frame, we only dumpy the samples between the triggers. */
              nconsumed += next_trigger;
              // Adjust phase offset, because we dumped samples.
              rotate_phase(&fine_freq_off[nconsumed], next_trigger);
            } else {
              /* There is no trigger at the beginning of this symbol, so we are
               * not at the beginning of the frame (where we would have to dump
               * the first 2 sync symbols). */
              // Copy rest of the data up to the trigger.
              for (int j = 0; j < next_trigger; ++j) { //TODO handle case: symbol_len>=next_trigger>fft_len
                // Calculate the modulated signal for fine frequency correction.
                gr_complex fine_freq_corr = std::polar((float)1.0, d_phase);
                for (int n = 0; n < d_n; ++n) {
                  ((gr_complex *) output_items[n])[nwritten + j] =
                          ((const gr_complex *) input_items[3 + n])[nconsumed + j] * fine_freq_corr;
                }
                // Rotate phase for fine frequency correction.
                rotate_phase(&fine_freq_off[nconsumed], 1);
              }
              // Add info tag 'out_of_sync'.
              add_item_tag(0,
                           (nitems_written(0) + nwritten)/d_fft_len,
                           pmt::mp("lost_sync"),
                           pmt::from_long(0));
              // Fill the rest of the symbol with zeros.
              for (int j = next_trigger; j < d_fft_len; ++j) {
                for (int n = 0; n < d_n; ++n) {
                  ((gr_complex *) output_items[n])[nwritten + j] = 0.0;
                }
              }
              nconsumed += next_trigger;
              nwritten += d_fft_len;
            }
            break;
          }
        }
        if (next_trigger == 0){
          // No triggers were found in the next symbol, check if there is a trigger at the beginning.
          if (trigger[nconsumed] == 1){
            // There is a trigger at the beginning of the symbol.
            // We dump the first of 2 (sync) symbols.
            nconsumed += d_symbol_len;
            d_sync_sym_count = 1;
            // Adjust phase offset, because we dumped samples.
            rotate_phase(&fine_freq_off[nconsumed], d_symbol_len);
          } else if(d_sync_sym_count == 1){
            // There is no trigger at the beginning of the symbol but there was one at the beginning of the previous symbol.
            // We dump the second of 2 (sync) symbols.
            nconsumed += d_symbol_len;
            d_sync_sym_count = 2;
            // Adjust phase offset, because we dumped samples.
            rotate_phase(&fine_freq_off[nconsumed], d_symbol_len);
          } else {
            // Check if this is the first symbol after the 2 sync symbols.
            if (d_sync_sym_count == 2){
              // Set start tag on this symbol.
              add_item_tag(0,
                           (nitems_written(0) + nwritten)/d_fft_len,
                           pmt::mp("start_of_frame"),
                           pmt::from_long(0));
              // Reset sync symbol counter.
              d_sync_sym_count = 0;
            }
            // There is no trigger at the beginning of the current and the previous symbol.
            // We dump the cyclic prefix and copy the rest.
            nconsumed += d_cp_len;
            // Adjust phase offset, because we dumped samples.
            rotate_phase(&fine_freq_off[nconsumed], d_cp_len);
            for (int j = 0; j < d_fft_len; ++j) {
              // Calculate the modulated signal for fine frequency correction.
              gr_complex fine_freq_corr = std::polar((float)1.0, d_phase);
              for (int n = 0; n < d_n; ++n) {
                ((gr_complex *) output_items[n])[nwritten + j] =
                        ((const gr_complex *) input_items[2 + n])[nconsumed + j] * fine_freq_corr;
              }
              // Rotate phase for fine frequency correction.
              rotate_phase(&fine_freq_off[nconsumed], 1);
            }
            nconsumed += d_fft_len;
            nwritten += d_fft_len;
          }
        }
      }
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (nconsumed);

      // Tell runtime system how many output items we produced.
      return noutput_items; // noutput_items = nwritten/d_symbol_len
    }

  } /* namespace digital */
} /* namespace gr */


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

#include <boost/format.hpp>
#include <gnuradio/io_signature.h>
#include "mimo_ofdm_synchronizer_fbcvc_impl.h"

using namespace boost;

namespace gr {
  namespace digital {

    mimo_ofdm_synchronizer_fbcvc::sptr
    mimo_ofdm_synchronizer_fbcvc::make(uint16_t n,
                                       uint32_t fft_len,
                                       uint32_t cp_len,
                                       const std::vector<gr_complex> &sync_symbol1,
                                       const std::vector<gr_complex> &sync_symbol2)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_synchronizer_fbcvc_impl(n,
                                               fft_len,
                                               cp_len,
                                               sync_symbol1,
                                               sync_symbol2));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_synchronizer_fbcvc_impl::mimo_ofdm_synchronizer_fbcvc_impl(
            uint16_t n,
            uint32_t fft_len,
            uint32_t cp_len,
            const std::vector<gr_complex> &sync_symbol1,
            const std::vector<gr_complex> &sync_symbol2)
      : gr::block("mimo_ofdm_synchronizer_fbcvc",
              gr::io_signature::make3(3+n, 3+n, sizeof(float), sizeof(unsigned char), sizeof(gr_complex)),
              gr::io_signature::make(n, n, sizeof(gr_complex) * fft_len)),
        d_n(n),
        d_fft_len(fft_len),
        d_cp_len(cp_len),
        d_on_frame(false),
        d_first_data_symbol(false),
        d_phase(0.),
        d_corr_v(sync_symbol2),
        d_first_active_carrier(0),
        d_last_active_carrier(sync_symbol2.size()-1)
    {
      // Check if both sync symbols have the length fft_len
      if (sync_symbol1.size() != sync_symbol2.size()) {
        throw std::invalid_argument("sync symbols must have equal length.");
      }
      if (sync_symbol1.size() != fft_len) {
        throw std::invalid_argument("sync symbols must have length fft_len.");
      }
      d_symbol_len = d_fft_len + d_cp_len;
      // Set up coarse freq estimation info
      // Allow all possible values:
      d_max_neg_carr_offset = -d_first_active_carrier;
      d_max_pos_carr_offset = d_fft_len - d_last_active_carrier - 1;
      d_rec_sync_symbol1 = std::vector<gr_complex> (fft_len, 0.0);
      d_rec_sync_symbol2 = std::vector<gr_complex> (fft_len, 0.0);

      set_tag_propagation_policy(TPP_DONT);
      set_min_noutput_items(2);
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
        ninput_items_required[i] = (noutput_items)*d_symbol_len;
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
    mimo_ofdm_synchronizer_fbcvc_impl::get_carr_offset(const std::vector<gr_complex> &sync_sym1,
                                                       const std::vector<gr_complex> &sync_sym2)
    {
      // TODO FFT before estimation!!!
      int carr_offset = 0;
      // Use Schmidl & Cox method
      float Bg_max = 0;
      // g here is 2g in the paper
      for (int g = d_max_neg_carr_offset; g <= d_max_pos_carr_offset; g += 2) {
        gr_complex tmp = gr_complex(0, 0);
        for (unsigned int k = 0; k < d_fft_len; k++) {
          if (d_corr_v[k] != gr_complex(0, 0)) {
            tmp += std::conj(sync_sym1[k+g]) * std::conj(d_corr_v[k]) * sync_sym2[k+g];
          }
        }
        if (std::abs(tmp) > Bg_max) {
          Bg_max = std::abs(tmp);
          carr_offset = g;
        }
      }
      return carr_offset;
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

      uint32_t symbol_count = 0;

      while (symbol_count < noutput_items) {
        // Search for the last trigger in the current symbol.
        uint32_t trigger_pos = 0;
        bool found_trigger = false;
        for (unsigned int k = 0; k < d_symbol_len; ++k){
          if(trigger[nconsumed+k] == 1){
            trigger_pos = std::max(trigger_pos, k);
            found_trigger = true;
          }
        }
        // State machine.
        if(found_trigger){
          GR_LOG_INFO(d_logger, format("-Rel %d, Abs %d. Found trigger at %d")%nconsumed %(nconsumed+nitems_read(0)) %trigger_pos);
          // Check if we can copy one fft_len vector before the trigger arrives.
          if(d_on_frame && (trigger_pos >= d_fft_len)){
            GR_LOG_INFO(d_logger, format("--Found trigger. Copy symbol before trigger."));
            // Copy this symbol before we dump the samples until the next trigger.
            // TODO copy symbol with rotator and fine freq correction multiplication.
            for (unsigned int i = 0; i < d_fft_len; ++i) {
              for (int n = 0; n < d_n; ++n) {
                ((gr_complex *) output_items[n])[nwritten + i] =
                        ((const gr_complex *) input_items[3 + n])[nconsumed + i];
              }
            }
            if(d_first_data_symbol) {
              // Set start tag on this fft vector.
              for (int n = 0; n < d_n; ++n) {
                add_item_tag(n,
                             nitems_written(0) + nwritten / d_fft_len,
                             pmt::string_to_symbol("start"),
                             pmt::from_long(0));
              }
              d_first_data_symbol = false;
            }
            nwritten += d_fft_len;
          }
          // Dump all samples of the new symbol until the trigger pos.
          nconsumed += trigger_pos;
          d_on_frame = false;
          if(nconsumed == 0){
            GR_LOG_INFO(d_logger, format("--Found trigger. Nconsumed == 0."));
            // The input buffer contains both sync symbols. We can process them if they are not interrupted by new triggers.
            // Search for triggers in the range of the sync symbols.
            trigger_pos = 0;
            for (unsigned int k = 1; k < 2*d_symbol_len; ++k){
              if(trigger[nconsumed+k] == 1) {
                trigger_pos = std::max(trigger_pos, k);
              }
            }
            if(trigger_pos == 0){
              GR_LOG_INFO(d_logger, format("---The sync syms are free of triggers."));
              // The sync syms dont get interrupted. Process them.
              // TODO Coarse freq estimation.
              nconsumed += 2*d_symbol_len;
              d_on_frame = true;
              d_first_data_symbol = true;
              // Consumed 2 symbols and adjust loop counter manually.
              symbol_count += 2;
              continue;
            } else {
              GR_LOG_INFO(d_logger, format("---The sync syms get interrupted."));
              // The sync syms get interrupted.
              // Start new synced buffer, which has min length of 2 sync symbols.
              break;
            }
          } else{
            GR_LOG_INFO(d_logger, format("--Found trigger. Nconsumed =%d.")%nconsumed);
            // Start new synced buffer, which has min length of 2 sync symbols.
            break;
          }
        } else {
          GR_LOG_INFO(d_logger, format("-Rel %d, abs %d. No trigger found.")%nconsumed %(nconsumed+nitems_read(0)));
          // We did not find a trigger.
          if(d_on_frame){
            GR_LOG_INFO(d_logger, format("-- On frame."));
            // We were on a frame and didn't find a trigger in this symbol.
            // Copy symbol and dump cyclic_prefix.

            // TODO copy symbol with rotator and fine freq correction multiplication.
            for (unsigned int i = 0; i < d_fft_len; ++i) {
              for (int n = 0; n < d_n; ++n) {
                ((gr_complex *) output_items[n])[nwritten + i] =
                        ((const gr_complex *) input_items[3 + n])[nconsumed + i];
              }
            }

            if(d_first_data_symbol) {
              // Set start tag on this fft vector.
              for (int n = 0; n < d_n; ++n) {
                add_item_tag(n,
                             nitems_written(0) + nwritten / d_fft_len,
                             pmt::string_to_symbol("start"),
                             pmt::from_long(0));
              }
              d_first_data_symbol = false;
            }
            nconsumed += d_symbol_len;
            nwritten += d_fft_len;
          } else {
            GR_LOG_INFO(d_logger, format("-- Not on frame."));
            // We were not on a frame and didn't find a trigger.
            // Dump whole symbol.
            nconsumed += d_symbol_len;
          }
        }
        symbol_count++;
      }

      consume_each (nconsumed);
      return nwritten/d_fft_len;
    }

  } /* namespace digital */
} /* namespace gr */


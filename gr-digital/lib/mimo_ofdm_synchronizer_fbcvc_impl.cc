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
#include <numeric>
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
                                       const std::vector<gr_complex> &sync_symbol2,
                                       const std::string &start_key)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_synchronizer_fbcvc_impl(n,
                                               fft_len,
                                               cp_len,
                                               sync_symbol1,
                                               sync_symbol2,
                                               start_key));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_synchronizer_fbcvc_impl::mimo_ofdm_synchronizer_fbcvc_impl(
            uint16_t n,
            uint32_t fft_len,
            uint32_t cp_len,
            const std::vector<gr_complex> &sync_symbol1,
            const std::vector<gr_complex> &sync_symbol2,
            const std::string &start_key)
      : gr::block("mimo_ofdm_synchronizer_fbcvc",
              gr::io_signature::make3(3+n, 3+n, sizeof(float), sizeof(unsigned char), sizeof(gr_complex)),
              gr::io_signature::make(n, n, sizeof(gr_complex) * fft_len)),
        d_n(n),
        d_fft_len(fft_len),
        d_cp_len(cp_len),
        d_symbol_len(fft_len+cp_len),
        d_on_frame(false),
        d_sync_read(false),
        d_first_data_symbol(false),
        d_phase(0.),
        d_carrier_freq_offset(0),
        d_first_active_carrier(0),
        d_last_active_carrier(sync_symbol2.size()-1),
        d_corr_v(sync_symbol2),
        d_start_key(pmt::string_to_symbol(start_key))
    {
      // Check if both sync symbols have the length fft_len
      if (sync_symbol1.size() != sync_symbol2.size()) {
        throw std::invalid_argument("Sync symbols must have equal length.");
      }
      if (sync_symbol1.size() != fft_len) {
        throw std::invalid_argument("Sync symbols must have length fft_len.");
      }

      // Set index of first and last active carrier.
      for (unsigned int i = 0; i < d_fft_len; i++) {
        if (sync_symbol1[i] != gr_complex(0, 0)) {
          d_first_active_carrier = i;
          break;
        }
      }
      for (int i = d_fft_len-1; i >= 0; i--) {
        if (sync_symbol1[i] != gr_complex(0, 0)) {
          d_last_active_carrier = i;
          break;
        }
      }

      /* Allocate buffer for FFT calculations.
       * (Required for integer carrier frequency offset estimation) */
      d_fft = new fft::fft_complex(d_fft_len, true, 1);
      d_rec_sync_symbol1 = std::vector<gr_complex> (fft_len, 0.0);
      d_rec_sync_symbol2 = std::vector<gr_complex> (fft_len, 0.0);

      // Set up coarse freq estimation parameters.
      d_max_neg_carr_offset = -d_first_active_carrier;
      d_max_pos_carr_offset = d_fft_len - d_last_active_carrier - 1;
      // Carrier offsets must be even.
      if (d_max_neg_carr_offset % 2)
        d_max_neg_carr_offset++;
      if (d_max_pos_carr_offset % 2)
        d_max_pos_carr_offset--;

      // Calculate differential modulated PN sequence
      for (unsigned int i = 0; i < d_fft_len; i++) {
        if (sync_symbol1[i] == gr_complex(0, 0)) {
          d_corr_v[i] = gr_complex(0, 0);
        } else {
          d_corr_v[i] /= sync_symbol1[i];
        }
      }

      // Set block behavior.
      set_tag_propagation_policy(TPP_DONT);
      // Ensure that a minimum of 2 OFDM symbols fit into the input buffer.
      set_output_multiple(2);
    }

    /*
     * Our virtual destructor.
     */
    mimo_ofdm_synchronizer_fbcvc_impl::~mimo_ofdm_synchronizer_fbcvc_impl()
    {
      // Delete allocation for FFT.
      delete d_fft;
    }

    void
    mimo_ofdm_synchronizer_fbcvc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      for (int i = 0; i < 3+d_n; ++i) {
        ninput_items_required[i] = (noutput_items)*d_symbol_len;
      }
    }

    uint32_t
    mimo_ofdm_synchronizer_fbcvc_impl::find_trigger(const unsigned char *trigger,
                                                    uint32_t start,
                                                    uint32_t end) {
      uint32_t trigger_pos = end;
      for (unsigned int k = start; k < end; ++k){ //TODO better search possible?
        if(trigger[k] != 0){
          trigger_pos = k;
          break;
        }
      }
      return trigger_pos;
    }

    void
    mimo_ofdm_synchronizer_fbcvc_impl::rotate_phase(const float *fine_freq_off,
                                                    uint16_t rotation_length) {
      d_phase += 2.0*std::accumulate(fine_freq_off, &fine_freq_off[rotation_length], 0.0)/d_fft_len;
      // Place the phase in [0, pi).
      d_phase = std::fmod(d_phase, 2.*M_PI);
    }

    void
    mimo_ofdm_synchronizer_fbcvc_impl::extract_symbols(gr_vector_const_void_star &input_items,
                                                       uint32_t input_offset,
                                                       gr_vector_void_star &output_items,
                                                       uint32_t output_offset,
                                                       uint32_t num_symbols) {
      const float *fine_freq_off = (const float *) input_items[0];

      for (unsigned int j = 0; j < num_symbols; ++j) {
        // Copy symbol.
        for (unsigned int i = 0; i < d_fft_len; ++i) {
          for (int n = 0; n < d_n; ++n) {
            ((gr_complex *) output_items[n])[output_offset + j * d_fft_len + i] =
                    ((const gr_complex *) input_items[3 + n])[input_offset + j * d_symbol_len + i] *
                    std::polar((float) 1.0, -d_phase);
          }
          // Rotate phase for fractional frequency correction.
          rotate_phase(&fine_freq_off[input_offset + j * d_symbol_len + i], 1);
        }
        // Rotate phase over cp.
        rotate_phase(&fine_freq_off[input_offset + j * d_symbol_len + d_fft_len], d_cp_len);
      }
    }

    int
    mimo_ofdm_synchronizer_fbcvc_impl::get_carr_offset(gr_vector_const_void_star &input_items,
                                                       uint32_t input_offset)
    {
      const gr_complex *ref_sig = (const gr_complex *) input_items[2];
      const float *fine_freq_off = (const float *) input_items[0];
      // Calculate FFT of the (fine frequency corrected) sync symbol 1.
      for (unsigned int i = 0; i < d_fft_len; ++i) {
        d_fft->get_inbuf()[i] = ref_sig[i]*std::polar((float)1.0, -d_phase);
        // Rotate phase for fractional frequency correction.
        rotate_phase(&fine_freq_off[input_offset+i], 1);
      }
      d_fft->execute();
      // Save FFT vector to array.
      gr_complex sync_sym1_fft[d_fft_len];
      memcpy(sync_sym1_fft, d_fft->get_outbuf(), d_fft_len*sizeof(gr_complex));

      // Rotate CP between the sync symbols.
      rotate_phase(&fine_freq_off[input_offset + d_fft_len], d_cp_len); // TODO not necessary?

      // Calculate FFT of the (fine frequency corrected) sync symbol 2.
      for (unsigned int i = 0; i < d_fft_len; ++i) {
        d_fft->get_inbuf()[i] = ref_sig[input_offset+i+d_symbol_len]*std::polar((float)1.0, -d_phase);
        // Rotate phase for fractional frequency correction.
        rotate_phase(&fine_freq_off[input_offset+i+d_symbol_len], 1);
      }
      d_fft->execute();
      // Save FFT vector to array.
      gr_complex sync_sym2_fft[d_fft_len];
      memcpy(sync_sym2_fft, d_fft->get_outbuf(), d_fft_len*sizeof(gr_complex));

      // Rotate CP after the second sync symbols.
      rotate_phase(&fine_freq_off[input_offset + d_symbol_len + d_fft_len], d_cp_len); // TODO not necessary?

      int carr_offset = 0;
      // Use method of Schmidl and Cox: "Robust Frequency and Timing Synchronization for OFDM."
      float Bg_max = 0;
      // g here is 2g in the paper.
      for (int g = d_max_neg_carr_offset; g <= d_max_pos_carr_offset; g += 2) {
        gr_complex tmp = gr_complex(0, 0);
        for (unsigned int k = 0; k < d_fft_len; k++) {
          if (d_corr_v[k] != gr_complex(0, 0)) {
            tmp += std::conj(sync_sym1_fft[(k+g+d_fft_len/2)%d_fft_len]) *
                    std::conj(d_corr_v[k]) * sync_sym2_fft[(k+g+d_fft_len/2)%d_fft_len];
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
      const unsigned char *trigger = (const unsigned char *) input_items[1];
      // Init state parameters.
      uint32_t noutput_samples = noutput_items*d_symbol_len;
      uint32_t nconsumed = 0;
      uint32_t nwritten = 0;

      if (d_on_frame){
        // We are currently on a frame.
        if (d_sync_read){
          // Search for new triggers in the buffer (indicating the beginning of a new frame).
          uint32_t trigger_pos = find_trigger(trigger, 0, noutput_samples);
          // Process symbols of the current frame.
          uint32_t num_syms = (trigger_pos+d_cp_len)/d_symbol_len; // The flooring is included in the implicit integer cast.
          // Copy symbols (without cp) to output buffer.
          extract_symbols(input_items, 0, output_items, 0, num_syms);

          // Set start tag if it is the first data symbol.
          if(d_first_data_symbol) {
            for (int n = 0; n < d_n; ++n) {
              add_item_tag(n, nitems_written(0) + 0, d_start_key, pmt::from_long(trigger_pos)); //TODO what to set as value?
              add_item_tag(n, nitems_written(0) + 0, pmt::string_to_symbol("carrier_freq_offset"),
                             pmt::from_long(d_carrier_freq_offset));
            }
            d_first_data_symbol = false;
          }
          nwritten = num_syms;
          // Check if we arrived at the end of the frame.
          if (trigger_pos < noutput_samples){
            // Detected start of new symbol.
            nconsumed = trigger_pos;
            // Fpr the new symbol, we have to read the sync symbols at first.
            d_sync_read = false;
          } else {
            // Frame continues with the next buffer.
            nconsumed = num_syms*d_symbol_len;
          }
        } else {
          // We are at the beginning of the frame (--> 2 sync symbols).
          // Read sync symbols and estimate integer carrier frequency offset.
          //d_phase = 0.0;
          d_carrier_freq_offset = get_carr_offset(input_items, 0);
          nconsumed = 2*d_symbol_len;
          d_sync_read = true;
          d_first_data_symbol = true;
        }
      } else {
        // We are not on a frame.
        // Search for a frame in the input stream.
        uint32_t trigger_pos = find_trigger(trigger, 0, noutput_samples);
        if (trigger_pos < noutput_samples){
          // Found a trigger. Dump samples until here.
          nconsumed = trigger_pos;
          // Now we are on a frame.
          d_on_frame = true;
          d_sync_read = false;
        } else {
          // Found no trigger.
          nconsumed = noutput_samples;
        }
      }
      consume_each(nconsumed);
      return nwritten;
    }

  } /* namespace digital */
} /* namespace gr */


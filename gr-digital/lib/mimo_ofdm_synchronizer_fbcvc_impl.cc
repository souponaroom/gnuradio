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
        d_packet_id(0),
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
      //set_min_noutput_items(2);
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
      for (int i = 0; i < rotation_length; ++i) { // TODO replace for loop with multiplication
        d_phase += fine_freq_off[i] * 2.0 / d_fft_len;
      }
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
        // Rotate phase over cp.
        rotate_phase(&fine_freq_off[input_offset + j * d_symbol_len], d_cp_len);
        // Copy symbol.
        for (unsigned int i = 0; i < d_fft_len; ++i) {
          for (int n = 0; n < d_n; ++n) {
            ((gr_complex *) output_items[n])[output_offset + j * d_fft_len + i] =
                    ((const gr_complex *) input_items[3 + n])[input_offset + j * d_symbol_len + d_cp_len + i];
                    //* std::polar((float) 1.0, -d_phase); TODO enable fine freq correction
          }
          // Rotate phase for fractional frequency correction.
          rotate_phase(&fine_freq_off[input_offset + j * d_symbol_len + d_cp_len * i], 1);
        }
      }
    }

    int
    mimo_ofdm_synchronizer_fbcvc_impl::get_carr_offset(const gr_complex *sync_sym1_fft,
                                                       const gr_complex *sync_sym2_fft)
    {
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
      const float *fine_freq_off = (const float *) input_items[0];
      const unsigned char *trigger = (const unsigned char *) input_items[1];
      const gr_complex *ref_sig = (const gr_complex *) input_items[2];
      uint32_t noutput_samples = noutput_items*d_symbol_len;
      uint32_t nconsumed = 0;
      uint32_t nwritten = 0;

      // Trigger points to beginning of cp!


      GR_LOG_INFO(d_logger, format("Noutput Items = %d")%noutput_items);
      if (d_on_frame){
        GR_LOG_INFO(d_logger, format("-- On frame."));
        if (d_sync_read){
              GR_LOG_INFO(d_logger, format("---- Already read sync."));
          // We already read the sync.
          // Search for triggers in the buffer.
          uint32_t trigger_pos = find_trigger(trigger, 0, noutput_samples);
          GR_LOG_INFO(d_logger, format("#######Trigger Pos %d")%trigger_pos);
          //trigger_pos = 160; // TODO Remove
          // Process symbols of the current frame.
          uint32_t num_syms = (trigger_pos+d_cp_len)/d_symbol_len; // The flooring is included in the implicit integer cast.
          //todo check if symbol_len*num_syms is in input buffer
          //extract_symbols(input_items, 0, output_items, 0, num_syms); TODO enable
          // TODO disable this
          const gr_complex *in1 = (const gr_complex *) input_items[3];
          const gr_complex *in2 = (const gr_complex *) input_items[4];
          gr_complex *out1 = (gr_complex*)output_items[0];
          gr_complex *out2 = (gr_complex*)output_items[1];
          GR_LOG_INFO(d_logger, format("---- Write %d syms at %d (packet %d)")%num_syms %(nitems_written(0)) %(d_packet_id));
          for (int i = 0; i < num_syms; ++i) {
            GR_LOG_INFO(d_logger, format("------ Read at %d")%(1.0*(nitems_read(0)+i*d_symbol_len)/d_symbol_len));
            memcpy(&out1[i*d_fft_len], &in1[d_cp_len+i*d_symbol_len], sizeof(gr_complex)*d_fft_len);
            memcpy(&out2[i*d_fft_len], &in2[d_cp_len+i*d_symbol_len], sizeof(gr_complex)*d_fft_len);
          }

          // Set start tag if it is the first data symbol.
          if(d_first_data_symbol) {
            for (int n = 0; n < d_n; ++n) {
              add_item_tag(n, nitems_written(0) + 0, d_start_key, pmt::from_long(noutput_samples));
            }
            d_first_data_symbol = false;
          }
          nwritten = num_syms;
          // Check if we were interrupted by a trigger.
          if (1){
          //if (trigger_pos < noutput_samples){
            GR_LOG_INFO(d_logger, format("------ Found new trigger in buffer at %d. (noutpuots samples %d)")%(trigger_pos) %noutput_samples);
            // Found a trigger.
            nconsumed = 160;//trigger_pos;
            d_sync_read = false;
            d_packet_id++;
          } else {
            GR_LOG_INFO(d_logger, format("------ No trigger in buffer."));
            // No triggers in this buffer. Call new work() to get new symbols of the current frame.
            nconsumed = std::min(num_syms*d_symbol_len, (uint32_t)noutput_samples); //todo necessary?
          }
        } else {
          // We are at the beginning of the frame.
          GR_LOG_INFO(d_logger, format("---- Read sync syms"));
          // Read sync symbols first.
          // TODO est freq (carr+int) offset by methods
          nconsumed = 2*d_symbol_len;
          d_sync_read = true;
          d_first_data_symbol = true;
        }
      } else {
        GR_LOG_INFO(d_logger, format("-- Not on frame."));
        // We are not on a frame.
        // Search for a trigger in the input stream.
        uint32_t trigger_pos = find_trigger(trigger, 0, noutput_samples);
        if (trigger_pos < noutput_samples){
          GR_LOG_INFO(d_logger, format("---- Found a trigger at pos %d.")%(nitems_read(0)+trigger_pos));
          // Found a trigger. Dump samples until here.
          nconsumed = trigger_pos;
          // Now we are on a frame.
          d_on_frame = true;
          d_sync_read = false;
        } else {
          GR_LOG_INFO(d_logger, format("---- Found no trigger."));
          // Found no trigger.
          nconsumed = noutput_samples;
        }
      }
      GR_LOG_INFO(d_logger, format("Nconsumed %d, Nwritten %d")%nconsumed %nwritten);

























      // OLD APPROACH

//      uint32_t symbol_count = 0;
//
//      while (symbol_count < noutput_items) {
//        // Search for the last trigger (indicating the detection of a new frame) in the current symbol.
//        uint32_t trigger_pos = 0;
//        bool found_trigger = false;
//        for (unsigned int k = d_symbol_len; k >= 1; --k){ //TODO better search possible?
//          if(trigger[nconsumed+k-1] == 1){
//            trigger_pos = k-1;
//            found_trigger = true;
//            break;
//          }
//        }
//        // State machine.
//        if(found_trigger){
//          // We detected the start of a new frame in the current symbol.
//          // Check if we can copy one vector (with fft_len elements) before the trigger arrives.
//          if(d_on_frame && (trigger_pos >= d_fft_len)){
//            // Copy this symbol before we dump the samples until the next trigger.
//            for (unsigned int i = 0; i < d_fft_len; ++i) {
//              // Copy symbol.
//              for (int n = 0; n < d_n; ++n) {
//                ((gr_complex *) output_items[n])[nwritten + i] = ((const gr_complex *) input_items[3 + n])[nconsumed + i]*std::polar((float)1.0, -d_phase);
//              }
//              // Rotate phase for fractional frequency correction.
//              rotate_phase(&fine_freq_off[nconsumed+i], 1);
//            }
//            // Check if this symbol was the first data symbol (after sync syms) of the frame.
//            if(d_first_data_symbol) {
//              // Set 'start' tag on this FFT vector.
//              for (int n = 0; n < d_n; ++n) {
//                add_item_tag(n,
//                             nitems_written(0) + nwritten / d_fft_len,
//                             d_start_key,
//                             pmt::from_long(0));
//              }
//              d_first_data_symbol = false;
//            }
//            nwritten += d_fft_len;
//            // Dump rest of the samples until the trigger pos and rotate phase.
//            rotate_phase(&fine_freq_off[nconsumed+d_fft_len], trigger_pos - d_fft_len);
//          } else {
//            // Dump samples before the next trigger.
//            rotate_phase(&fine_freq_off[nconsumed], trigger_pos);
//          }
//          // Dump all samples until the trigger.
//          nconsumed += trigger_pos;
//          d_on_frame = false;
//          // We are finished with the symbol before the trigger.
//
//          /* The first two symbols after the trigger should be sync symbols.
//           * Check if the next two symbol vectors are completely in the
//           * current buffer. */
//          if(nconsumed == 0){
//            /* The input buffer contains both sync symbols.
//             * We can process them if they are not interrupted by any new triggers. */
//            // Search for triggers in the range of the sync symbols.
//            trigger_pos = 0;
//            for (unsigned int k = 2*d_symbol_len; k >= 2; --k){
//              if(trigger[k-1] == 1) {
//                trigger_pos = k-1;
//                break;
//              }
//            }
//            if(trigger_pos == 0){
//              // The sync syms dont get interrupted. Process them.
//              // Estimate integer carrier frequency offset. TODO integrate this stuff into the method get_carr_offset
//
//              // Calculate FFT of the (fine frequency corrected) sync symbol 1.
//              for (unsigned int i = 0; i < d_fft_len; ++i) {
//                d_fft->get_inbuf()[i] = ref_sig[i]*std::polar((float)1.0, -d_phase);
//                // Rotate phase for fractional frequency correction.
//                rotate_phase(&fine_freq_off[i], 1);
//              }
//              d_fft->execute();
//              // Save FFT vector to array.
//              gr_complex sync_sym1_fft[d_fft_len];
//              memcpy(sync_sym1_fft, d_fft->get_outbuf(), d_fft_len*sizeof(gr_complex));
//
//              // Rotate CP between the sync symbols.
//              rotate_phase(&fine_freq_off[d_fft_len], d_cp_len); // TODO not necessary?
//
//              // Calculate FFT of the (fine frequency corrected) sync symbol 2.
//              for (unsigned int i = 0; i < d_fft_len; ++i) {
//                d_fft->get_inbuf()[i] = ref_sig[i+d_symbol_len]*std::polar((float)1.0, -d_phase);
//                // Rotate phase for fractional frequency correction.
//                rotate_phase(&fine_freq_off[i+d_symbol_len], 1);
//              }
//              d_fft->execute();
//              // Save FFT vector to array.
//              gr_complex sync_sym2_fft[d_fft_len];
//              memcpy(sync_sym2_fft, d_fft->get_outbuf(), d_fft_len*sizeof(gr_complex));
//
//              // Estimate carrier frequency offset.
//              d_carrier_freq_offset = get_carr_offset(sync_sym1_fft, sync_sym2_fft);
//              // Rotate CP after the second sync symbol.
//              rotate_phase(&fine_freq_off[d_symbol_len+d_fft_len], d_cp_len);
//              nconsumed = 2*d_symbol_len;
//              d_on_frame = true;
//              d_first_data_symbol = true;
//              // Consumed 2 symbols and adjust loop counter manually.
//              symbol_count += 2;
//              continue;
//            } else {
//              GR_LOG_INFO(d_logger, format("The sync symbols get interrupted by a new trigger signal."));
//              /* The sync syms get interrupted. The frame detection of the
//               * timing synchronization may too sensible. */
//              // Start new synced buffer, which has min length of 2 sync symbols.
//              break;
//            }
//          } else{
//            // Start new synced buffer, which has min length of 2 sync symbols.
//            break;
//          }
//        } else {
//          // We did not find a trigger.
//          if(d_on_frame){
//            // We were on a frame and didn't find a trigger in this symbol.
//            // Copy symbol and dump cyclic_prefix.
//            for (unsigned int i = 0; i < d_fft_len; ++i) {
//              for (int n = 0; n < d_n; ++n) {
//                ((gr_complex *) output_items[n])[nwritten + i] = ((const gr_complex *) input_items[3 + n])[nconsumed + i]*std::polar((float)1.0, -d_phase);
//              }
//              // Rotate phase for fractional frequency correction.
//              rotate_phase(&fine_freq_off[nconsumed+i], 1);
//            }
//            // Check if this symbol was the first data symbol (after sync syms) of the frame.
//            if(d_first_data_symbol) {
//              // Set start tag on this fft vector.
//              for (int n = 0; n < d_n; ++n) {
//                add_item_tag(n,
//                             nitems_written(0) + nwritten / d_fft_len,
//                             pmt::string_to_symbol("start"),
//                             pmt::from_long(0));
//                add_item_tag(n,
//                             nitems_written(0) + nwritten / d_fft_len,
//                             pmt::string_to_symbol("carrier_freq_offset"),
//                             pmt::from_long(d_carrier_freq_offset));
//              }
//              d_first_data_symbol = false;
//            }
//            // Rotate phase for fractional frequency correction.
//            rotate_phase(&fine_freq_off[nconsumed + d_fft_len], d_cp_len);
//            nconsumed += d_symbol_len;
//            nwritten += d_fft_len;
//          } else {
//            // We were not on a frame and didn't find a trigger.
//            // Dump whole symbol.
//            rotate_phase(&fine_freq_off[nconsumed], d_symbol_len);
//            nconsumed += d_symbol_len;
//          }
//        }
//        symbol_count++;
//      }



      consume_each (nconsumed);
      return nwritten;
    }

  } /* namespace digital */
} /* namespace gr */


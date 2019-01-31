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
#include "mimo_ofdm_channel_estimator_vcvc_impl.h"
#include <pmt/pmt.h>

using namespace boost;

namespace gr {
  namespace digital {

    mimo_ofdm_channel_estimator_vcvc::sptr
    mimo_ofdm_channel_estimator_vcvc::make(uint16_t m, uint16_t n,
                                           uint32_t fft_len,
                                           std::vector<std::vector<gr_complex> > pilot_symbols,
                                           std::vector<int> pilot_carriers,
                                           std::vector<int> occupied_carriers,
                                           const std::string &csi_key,
                                           const std::string &start_key)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_channel_estimator_vcvc_impl(m, n,
                                                   fft_len,
                                                   pilot_symbols,
                                                   pilot_carriers,
                                                   occupied_carriers,
                                                   csi_key,
                                                   start_key));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_channel_estimator_vcvc_impl::mimo_ofdm_channel_estimator_vcvc_impl(
            uint16_t m, uint16_t n,
            uint32_t fft_len, 
            std::vector<std::vector<gr_complex> > pilot_symbols, 
            std::vector<int> pilot_carriers,
            std::vector<int> occupied_carriers,
            const std::string &csi_key,
            const std::string &start_key)
      : gr::block("mimo_ofdm_channel_estimator_vcvc",
              gr::io_signature::make(n, n, sizeof(gr_complex)*fft_len),
              gr::io_signature::make(n, n, sizeof(gr_complex)*occupied_carriers.size())),
        d_m(m), d_n(n),
        d_fft_len(fft_len),
        d_fft_shift(fft_len/2),
        d_pilot_symbols(pilot_symbols),
        d_pilot_carriers(pilot_carriers),
        d_occupied_carriers(occupied_carriers),
        d_output_vlen(occupied_carriers.size()),
        d_csi_key(pmt::string_to_symbol(csi_key)),
        d_start_key(pmt::string_to_symbol(start_key)),
        d_last_tag_offset(0)
    {
      // Init CSI vector.
      d_channel_state = std::vector<std::vector<std::vector<gr_complex> > >
              (d_fft_len, std::vector<std::vector<gr_complex> > (n, std::vector<gr_complex> (m, 1.0)));
      // Monitor that FFT length is even.
      if (fft_len%2){
        throw std::invalid_argument((boost::format("FFT length %d must be even.") %fft_len).str());
      }
      // Set tag propagation policy.
      set_tag_propagation_policy(TPP_ONE_TO_ONE);
    }

    /*
     * Our virtual destructor.
     */
    mimo_ofdm_channel_estimator_vcvc_impl::~mimo_ofdm_channel_estimator_vcvc_impl()
    {
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      for (int i = 0; i < d_m; ++i) {
        ninput_items_required[i] = noutput_items+d_m-1;
      }
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::extract_payload_carriers(
            gr_vector_const_void_star &input_items,
            gr_vector_void_star &output_items,
            uint32_t length) {
      // Iterate over RX branches.
      for (int i = 0; i < d_n; ++i) {
        const gr_complex *in = (const gr_complex *) input_items[i];
        gr_complex *out = (gr_complex *) output_items[i];
        // Extract data carriers to output buffer and dump zero-carriers and pilot-carriers.
        for (unsigned int j = 0; j < length; ++j) {
          for (unsigned int k = 0; k < d_output_vlen; ++k) {
            out[j*d_output_vlen + k] = in[j*d_fft_len + d_occupied_carriers[k] + d_fft_shift];
          }
        }
      }
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::estimate_channel_state(
            gr_vector_const_void_star &input_items,
            uint32_t reading_offset,
            uint16_t correlation_offset) {
      // Iterate over all pilot carriers.
      for (unsigned int c = 0; c < d_pilot_carriers.size(); ++c) {
        // Iterate over N MIMO RX streams.
        for (int i = 0; i < d_n; ++i) {
          const gr_complex *in = (const gr_complex *) input_items[i];
          for (int j = 0; j < d_m; ++j) { // Iterate over MIMO pilot sequences of length M.
            /* Correlate received pilot sequence with transmitted sequence.
             * The result is the path coefficient h_ij for the current sub-carrier
             * and the current MIMO branch. */
            d_channel_state[d_pilot_carriers[c]+d_fft_shift][i][j]
                    = correlate_pilots(&in[reading_offset*d_fft_len + d_pilot_carriers[c] + d_fft_shift],
                    d_pilot_symbols[j],
                    d_fft_len,
                    correlation_offset);
          }
        }
      }
    }

    gr_complex
    mimo_ofdm_channel_estimator_vcvc_impl::correlate_pilots(const gr_complex *in,
                                                            std::vector<gr_complex> pilot,
                                                            uint32_t distance,
                                                            uint16_t pilot_offset) {
      gr_complex correlation = 0.0;
      gr_complex energy = 0.0;
      // Correlate the received pilot sequence with the transmitted one.
      for (unsigned int i = 0; i < d_m; ++i) {
        correlation += in[i*distance] * std::conj(pilot[(i+pilot_offset)%d_m]);
        energy += in[i*distance]*std::conj(in[i*distance]);
      }
      // Return normalized (assuming a normalized pilot sequence) correlation result.
      return correlation/(gr_complex)d_m;
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::interpolate_channel_state() {
      // Set the CSI for the edge carriers to the same as their nearest estimation.
      for (unsigned int c = 0; c < d_pilot_carriers[0]+d_fft_shift; ++c) {
        d_channel_state[c] = d_channel_state[d_pilot_carriers[0]+d_fft_shift];
      }
      for (unsigned int c = d_pilot_carriers[d_pilot_carriers.size()-1]+d_fft_shift+1; c < d_fft_len; ++c) {
        d_channel_state[c] = d_channel_state[d_pilot_carriers[d_pilot_carriers.size()-1]+d_fft_shift];
      }
      // Linear interpolation over the OFDM sub-carriers of the complex channel coefficients.
      for (unsigned int k = 0; k < d_pilot_carriers.size()-1; ++k) {
        for (unsigned int c = d_pilot_carriers[k]+d_fft_shift; c < d_pilot_carriers[k+1]+d_fft_shift; ++c) {
          for (int i = 0; i < d_n; ++i) {
            for (int j = 0; j < d_m; ++j) {
              d_channel_state[c][i][j] =
                      d_channel_state[d_pilot_carriers[k]+d_fft_shift][i][j] +
                      ((gr_complex)(c-d_pilot_carriers[k]-d_fft_shift)/(gr_complex)(d_pilot_carriers[k+1]-d_pilot_carriers[k]))*
                              (d_channel_state[d_pilot_carriers[k+1]+d_fft_shift][i][j] -
                                      d_channel_state[d_pilot_carriers[k]+d_fft_shift][i][j]);
            }
          }
        }
      }
    }

    pmt::pmt_t
    mimo_ofdm_channel_estimator_vcvc_impl::generate_csi_pmt() {
      // Assign the channel state vector to a PMT vector. Only take the occupied carriers into account.
      pmt::pmt_t csi_pmt = pmt::make_vector(d_output_vlen, pmt::make_vector(d_m, pmt::make_c32vector(d_m, d_channel_state[0][0][0])));
      // Frequency dimension.
      for (unsigned int k = 0; k < d_output_vlen; ++k) {
        pmt::pmt_t csi_per_carrier = pmt::make_vector(d_n, pmt::make_c32vector(d_m, d_channel_state[0][0][0]));
        // RX space dimension.
        for (int i = 0; i < d_n; ++i){
          pmt::pmt_t csi_line_vector = pmt::make_c32vector(d_m, d_channel_state[0][0][0]);
          // TX space dimension.
          for (int j = 0; j < d_m; ++j) {
            pmt::c32vector_set(csi_line_vector, j, d_channel_state[d_occupied_carriers[k]+d_fft_len/2][i][j]);
          }
          pmt::vector_set(csi_per_carrier, i, csi_line_vector);
        }
        pmt::vector_set(csi_pmt, k, csi_per_carrier);
      }
      return csi_pmt;
    }

    int
    mimo_ofdm_channel_estimator_vcvc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      add_item_tag(0, nitems_written(0), pmt::string_to_symbol("est_buffer"), pmt::from_long(noutput_items));
      /* Copy occupied OFDM sub-carriers to output buffer.
       * (Neither the zero-carriers nor the pilot carriers) */
      extract_payload_carriers(input_items, output_items, noutput_items);

      // Read start tags to sync the pilot correlation.
      std::vector <gr::tag_t> start_tags;
      get_tags_in_window(start_tags, 0, 0, noutput_items+d_m-1, d_start_key);

      uint32_t tag_offset_correction = 0;
      if (nitems_read(0) >= d_last_tag_offset){
                tag_offset_correction = nitems_read(0)-d_last_tag_offset;
      }
      uint16_t tag_index = 0;

      // Iterate over OFDM symbols.
      for (int s = 0; s < noutput_items; ++s) {
        if(start_tags.size() > 0 && nitems_read(0)+s == start_tags[tag_index].offset){
          // Update current start offset.
          tag_offset_correction = d_m - (s%d_m);
          if(tag_index < start_tags.size()-1){
            tag_index++;
          }
        }
        // Experimental feature
        if(start_tags.size() > 0 && start_tags[tag_index].offset-(d_m-1) <= nitems_read(0)+s && nitems_read(0)+s < start_tags[tag_index].offset){
          // This is the last symbol of the frame.
          // Use old estimation.
        } else {
          // Estimate the complex channel coefficient of all pilot carriers (of all MIMO branches).
          estimate_channel_state(input_items, s, (s + tag_offset_correction) % d_m);
          /* We have estimated the CSI for the pilot carriers.
           * Now, lets interpolate over all remaining OFDM sub-carriers. */
          interpolate_channel_state();
        }
        /* Now we have individual CSI for each sub-carrier of each MIMO-branch.
         * Add tag with this CSI to the output vector.
         * All CSI for one time step is stored in one 3-dim vector which is tagged
         * to the output vector of the first MIMO branch. */
        add_item_tag(0, nitems_written(0) + s, d_csi_key, generate_csi_pmt());
      }

      // Update last_tag_offset.
      if (start_tags.size() > 0){
        d_last_tag_offset = start_tags[start_tags.size()-1].offset;
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


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
#include "mimo_ofdm_channel_estimator_vcvc_impl.h"
#include <pmt/pmt.h>

using namespace boost;

namespace gr {
  namespace digital {

    mimo_ofdm_channel_estimator_vcvc::sptr
    mimo_ofdm_channel_estimator_vcvc::make(uint16_t n,
                                           uint32_t fft_len,
                                           std::vector<std::vector<gr_complex> > pilot_symbols,
                                           std::vector<int> pilot_carriers)
    {
      return gnuradio::get_initial_sptr
        (new mimo_ofdm_channel_estimator_vcvc_impl(n, fft_len, pilot_symbols, pilot_carriers));
    }

    /*
     * The private constructor
     */
    mimo_ofdm_channel_estimator_vcvc_impl::mimo_ofdm_channel_estimator_vcvc_impl(
            uint16_t n, 
            uint32_t fft_len, 
            std::vector<std::vector<gr_complex> > pilot_symbols, 
            std::vector<int> pilot_carriers)
      : gr::block("mimo_ofdm_channel_estimator_vcvc",
              gr::io_signature::make(n, n, sizeof(gr_complex)*fft_len),
              gr::io_signature::make(n, n, sizeof(gr_complex)*fft_len)),
        d_n(n),
        d_fft_len(fft_len),
        d_pilot_symbols(pilot_symbols),
        d_pilot_carriers(pilot_carriers)
    {
      d_channel_state = std::vector<std::vector<std::vector<gr_complex> > >
              (d_fft_len, std::vector<std::vector<gr_complex> > (n, std::vector<gr_complex> (n, 1.0)));
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
      ninput_items_required[0] = noutput_items+d_n-1;
    }

    gr_complex
    mimo_ofdm_channel_estimator_vcvc_impl::correlate_pilots(const gr_complex *in,
                                                            std::vector<gr_complex> pilot,
                                                            uint32_t distance,
                                                            uint16_t pilot_offset) {
      gr_complex correlation = 0.0;
      for (unsigned int i = 0; i < d_n; ++i) {
        correlation += in[i*distance] * std::conj(pilot[(i+pilot_offset)%d_n]);
      }
      // Return correlation after normalization.
      return correlation/(gr_complex)d_n;
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::interpolate_channel_state() {
      // Set the CSI for the edge carriers to the same as their nearest estimation.
      for (unsigned int c = 0; c < d_pilot_carriers[0]+d_fft_len/2; ++c) {
        d_channel_state[c] = d_channel_state[d_pilot_carriers[0]+d_fft_len/2];
      }
      for (unsigned int c = d_pilot_carriers[d_pilot_carriers.size()-1]+d_fft_len/2; c < d_fft_len; ++c) {
        d_channel_state[c] = d_channel_state[d_pilot_carriers[d_pilot_carriers.size()-1]+d_fft_len/2];
      }
      // Linear interpolation over the channel estimations of the OFDM carriers.
      for (unsigned int k = 0; k < d_pilot_carriers.size()-1; ++k) {
        for (unsigned int c = d_pilot_carriers[k]+d_fft_len/2; c < d_pilot_carriers[k+1]+d_fft_len/2; ++c) {
          for (int i = 0; i < d_n; ++i) {
            for (int j = 0; j < d_n; ++j) {
              d_channel_state[c][i][j] =
                      d_channel_state[d_pilot_carriers[k]+d_fft_len/2][i][j] +
                      ((gr_complex)(c-d_pilot_carriers[k]-d_fft_len/2)/(gr_complex)(d_pilot_carriers[k+1]-d_pilot_carriers[k]))*
                              (d_channel_state[d_pilot_carriers[k+1]+d_fft_len/2][i][j] -
                                      d_channel_state[d_pilot_carriers[k]+d_fft_len/2][i][j]);
            }
          }
        }
      }
    }

    int
    mimo_ofdm_channel_estimator_vcvc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      // Copy data to output.
      for (int i = 0; i < d_n; ++i) {
        const gr_complex *in = (const gr_complex *) input_items[i];
        gr_complex *out = (gr_complex *) output_items[i];
        memcpy(out, in, sizeof(gr_complex)*d_fft_len*noutput_items);
      }

      for (int s = 0; s < noutput_items; ++s) { // Iterate over OFDM symbols.
        for (unsigned int c = 0; c < d_pilot_carriers.size(); ++c) { // Iterate over pilot carriers
          for (int i = 0; i < d_n; ++i) { // Iterate over N MIMO input streams.
            const gr_complex *in = (const gr_complex *) input_items[i];
            //GR_LOG_DEBUG(d_logger, format("sym %d, carrier %d, stream %i: %d")%s %d_pilot_carriers[c] %i %in[d_pilot_carriers[c] + d_fft_len / 2 + d_fft_len * s]);
            for (int j = 0; j < d_n; ++j) { // Iterate over N MIMO pilot sequences.
              // Correlate received pilot symbols with reference. The result is the path coefficient h_ij.
              d_channel_state[d_pilot_carriers[c] + d_fft_len / 2][i][j] = correlate_pilots(
                      &in[d_pilot_carriers[c] + d_fft_len / 2 + d_fft_len * s],
                      d_pilot_symbols[j],
                      d_fft_len,
                      s%d_n); //TODO This only applies for the synced case. Otherwise, read 'start' tag!
            }
          }
        }
        // We have estimated the CSI for the pilot carriers. Now, lets interpolate over all OFDM carriers.
        interpolate_channel_state();

        // Assign the channel state vector to a PMT vector.
        pmt::pmt_t csi_pmt = pmt::make_vector(d_fft_len, pmt::make_vector(d_n, pmt::make_c32vector(d_n, d_channel_state[0][0][0])));
        for (unsigned int k = 0; k < d_fft_len; ++k) {
          pmt::pmt_t csi_per_carrier = pmt::make_vector(d_n, pmt::make_c32vector(d_n, d_channel_state[0][0][0]));
          for (int i = 0; i < d_n; ++i){
            pmt::pmt_t csi_line_vector = pmt::make_c32vector(d_n, d_channel_state[0][0][0]);
            for (int j = 0; j < d_n; ++j) {
              pmt::c32vector_set(csi_line_vector, j, d_channel_state[k][i][j]);
            }
            pmt::vector_set(csi_per_carrier, i, csi_line_vector);
          }
          pmt::vector_set(csi_pmt, k, csi_per_carrier);
        }
        // Add tag to output data.
        add_item_tag(0, nitems_written(0) + s, pmt::string_to_symbol(std::string("csi")), csi_pmt);
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


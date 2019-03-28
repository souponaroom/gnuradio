/* -*- c++ -*- */
/* 
 * Copyright 2019 Free Software Foundation, Inc.
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
#include "ofdm_snr_est_vcvc_impl.h"
#include <volk/volk.h>
#include <pmt/pmt.h>

#include <boost/format.hpp>
using namespace boost;

namespace gr {
  namespace digital {

    ofdm_snr_est_vcvc::sptr
    ofdm_snr_est_vcvc::make(uint16_t num_inputs, uint32_t fft_len,
                            std::vector<int> occupied_carriers,
                            std::vector<int> zero_carriers,
                            const std::string &snr_key,
                            uint16_t update_time, float averaging_length)
    {
      return gnuradio::get_initial_sptr
        (new ofdm_snr_est_vcvc_impl(num_inputs, fft_len, occupied_carriers,
                                    zero_carriers, snr_key,
                                    update_time, averaging_length));
    }

    /*
     * The private constructor
     */
    ofdm_snr_est_vcvc_impl::ofdm_snr_est_vcvc_impl(uint16_t num_inputs,
                                                   uint32_t fft_len,
                                                   std::vector<int> occupied_carriers,
                                                   std::vector<int> zero_carriers,
                                                   const std::string &snr_key,
                                                   uint16_t update_time,
                                                   float averaging_length)
      : gr::sync_block("ofdm_snr_est_vcvc",
              gr::io_signature::make(num_inputs, num_inputs, fft_len*sizeof(gr_complex)),
              gr::io_signature::make(num_inputs, num_inputs, fft_len*sizeof(gr_complex))),
        d_num_inputs(num_inputs), d_fft_len(fft_len),
        d_occupied_carriers(occupied_carriers), d_zero_carriers(zero_carriers),
        d_snr_key(pmt::string_to_symbol(snr_key)),
        d_update_time(update_time), d_averaging_length(averaging_length),
        d_snr(0.0)
    {
      d_snr = std::vector<float>(num_inputs, 2.0);
      unsigned int alignment = volk_get_alignment();
      d_mag_squared = (float *) volk_malloc(sizeof(float)*fft_len, alignment);
      set_output_multiple(update_time);
    }

    /*
     * Our virtual destructor.
     */
    ofdm_snr_est_vcvc_impl::~ofdm_snr_est_vcvc_impl()
    {
    }

    int
    ofdm_snr_est_vcvc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {

      for (int l = 0; l < noutput_items/d_update_time; ++l) {
        for (int i = 0; i < d_num_inputs; ++i) {

          const gr_complex *in = &((const gr_complex *) input_items[i])[l*d_fft_len*d_update_time];
          gr_complex *out = &((gr_complex *) output_items[i])[l*d_fft_len*d_update_time];

          // calc mag square of sub-carriers
          volk_32fc_magnitude_squared_32f(d_mag_squared, &in[i * d_fft_len], d_fft_len);
          // measure normalized energy of occupied sub-carriers
          float occupied_energy = 0;
          float noise_energy;
          for (int j = 0; j < d_occupied_carriers.size(); ++j) {
            occupied_energy += d_mag_squared[d_occupied_carriers[j]+d_fft_len/2];
          }
          occupied_energy /= d_occupied_carriers.size();
          for (int k = 0; k < d_zero_carriers.size(); ++k) {
            noise_energy += d_mag_squared[d_zero_carriers[k]+d_fft_len/2];
          }
          noise_energy /= d_zero_carriers.size();
          if(occupied_energy > noise_energy) {
            float current_snr = (occupied_energy - noise_energy) / noise_energy;
            d_snr[i] = (d_snr[i] * (d_averaging_length - 1) + current_snr) / d_averaging_length;
          } else {
            // This should not happen.
            GR_LOG_INFO(d_logger, format("Detected higher energy for zero carriers (%d) than"
                                                 "for occupied carriers (%d).")
                                  %noise_energy %occupied_energy);
          }
        }
        pmt::pmt_t snr_tag = pmt::make_f32vector(d_num_inputs, 0);
        for (int n = 0; n < d_num_inputs; ++n) {
          pmt::f32vector_set(snr_tag, n, d_snr[n]);
        }
        add_item_tag(0, nitems_written(0)+l*d_update_time, d_snr_key, snr_tag);
      }
      for (int n = 0; n < d_num_inputs; ++n) {
        const gr_complex *in = (const gr_complex *) input_items[n];
        gr_complex *out = (gr_complex *) output_items[n];
        // Copy input to output.
        memcpy(out, in, sizeof(gr_complex) * d_fft_len * noutput_items);
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace digital */
} /* namespace gr */


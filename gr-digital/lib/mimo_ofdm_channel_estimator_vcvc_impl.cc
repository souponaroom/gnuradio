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
const float  PI_F=3.14159265358979f;

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
        d_start_new_packet(true),
        d_time_delay(0),
        d_pilot_symbols(pilot_symbols),
        d_pilot_carriers(pilot_carriers),
        d_occupied_carriers(occupied_carriers),
        d_output_vlen(occupied_carriers.size()),
        d_csi_key(pmt::string_to_symbol(csi_key)),
        d_start_key(pmt::string_to_symbol(start_key)),
        d_correlation_offset(0)
    {
      // Init CSI vector.
      d_channel_state = std::vector<std::vector<std::vector<gr_complex> > >
              (d_fft_len, std::vector<std::vector<gr_complex> > (n, std::vector<gr_complex> (m, 1.0)));
      // Monitor that FFT length is even.
      if (fft_len%2){
        throw std::invalid_argument((boost::format("FFT length %d must be even.") %fft_len).str());
      }
      // Set tag propagation policy.
      set_tag_propagation_policy(TPP_DONT);
      
      set_history(d_n);
      set_output_multiple(d_n);
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
        for(int i=0; i<d_n; ++i)
        {
            ninput_items_required[i] = noutput_items;
        }
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::extract_payload_carriers(
            gr_vector_const_void_star &input_items,
            gr_vector_void_star &output_items,
            uint32_t offset,
            uint32_t length) {
      // Iterate over RX branches.
      for (int i = 0; i < d_n; ++i) {
        const gr_complex *in = (const gr_complex *) input_items[i];
        gr_complex *out = (gr_complex *) output_items[i];
        // Extract data carriers to output buffer and dump zero-carriers and pilot-carriers.
        for (unsigned int j = 0; j < length; ++j) {
          for (unsigned int k = 0; k < d_output_vlen; ++k) {
            out[j*d_output_vlen + k] = in[(offset+j)*d_fft_len + d_occupied_carriers[k] + d_fft_shift] *
                                       std::polar(1.f, (float)(d_occupied_carriers[k])*-2.f*PI_F*d_time_delay/d_fft_len);
          }
        }
      }
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::estimate_time_delay(
            gr_vector_const_void_star &input_items,
            uint32_t reading_offset) {
      gr_complex correlation;
      float max = 0;
      int max_index = -5;
      for (int p = -5; p<6; ++p)
      {
          correlation = 0;
          for(int q=0; q<d_n; ++q)//TODO fix this
          {
            const gr_complex *in = (const gr_complex *) input_items[q];
            for (unsigned int c = 0; c < d_pilot_carriers.size(); ++c) 
            {
                correlation += in[reading_offset*d_fft_len + d_pilot_carriers[c] + d_fft_shift] * 
                               std::polar(1.f, (float)(d_pilot_carriers[c])*-2.f*PI_F*p/d_fft_len); //TODO verallgemeinern!
            }
          }
          if(std::abs(correlation) > max) 
          {
            max = std::abs(correlation);
            max_index = p;
          }
      }
       //Iterate over all pilot carriers.
         //Iterate over N MIMO RX streams.
        for (int i = 0; i < d_n; ++i) {
          for (int j = 0; j < d_m; ++j) { // Iterate over MIMO pilot sequences of length M.

          }
        }
    }

    void
    mimo_ofdm_channel_estimator_vcvc_impl::estimate_time_delay2(
            gr_vector_const_void_star &input_items,
            uint32_t reading_offset) {
      float correlation;
      float min = d_n*d_pilot_carriers.size()*4;
      int min_index = 0;
      for (int p = 0; p<4; ++p)
      {
          correlation = 0;
          for(int q=0; q<d_n; ++q)//TODO fix this
          {
            const gr_complex *in = (const gr_complex *) input_items[q];
            for (unsigned int c = 1; c < d_pilot_carriers.size(); ++c) 
            {
                correlation += std::abs(std::arg(in[reading_offset*d_fft_len + d_pilot_carriers[c] + d_fft_shift] * 
                               std::polar(1.f, (float)(d_pilot_carriers[c])*2.f*PI_F*p/d_fft_len)) -  
                               std::arg(in[reading_offset*d_fft_len + d_pilot_carriers[c-1] + d_fft_shift] * 
                               std::polar(1.f, (float)(d_pilot_carriers[c-1])*2.f*PI_F*p/d_fft_len))); 
            }
          }
          //std::cout << "correlation " << p << ": " << correlation << std::endl;
          if(correlation < min) 
          {
            min = correlation;
            min_index = p;
          }
      }
      //std::cout << "delay2 " << min_index << std::endl;
      d_time_delay = min_index;
      std::cout << "time delay " << d_time_delay << std::endl;
      //d_time_delay = max_index;
       //Iterate over all pilot carriers.
         //Iterate over N MIMO RX streams.
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
                    correlation_offset) * 
                    std::polar(1.f, (float)(d_pilot_carriers[c])*-2.f*PI_F*d_time_delay/d_fft_len);
          }
        }
      }
            //std::cout << "Chanest estimation " << d_channel_state[d_pilot_carriers[0]+d_fft_shift][0][0] << std::endl;
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
        //std::cout << "chanest received pilots" << in[i*distance] << std::endl;
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

      // Search for next tag in inbuffer.
      std::vector <gr::tag_t> start_tags;
      get_tags_in_window(start_tags, 0, 0, noutput_items, d_start_key);
      uint32_t items_to_process;

      //std::cout<<"chanest, noutput_items " << noutput_items  << " tags " << start_tags.size() << ", new packet " << d_start_new_packet << std::endl;
      for(int i=0; i<start_tags.size();i++)
      {
       // std::cout << "tag position: " << start_tags[i].offset << std::endl;
      }
      if(d_start_new_packet)
      {
        //std::cout << "new packet " << std::endl;
        // We are at the beginning of a new packet. Sanity check that there is a tag at first position.
        if(start_tags.size() < 1 || start_tags[0].offset - nitems_read(0) != 0)
        {
            std::cout << "ERROR, no tag at beginning of packet." << std::endl;
        }
        items_to_process = d_n;
        /* Estimate the channel of the first N symbols of the new packet
         * without the use of history symbols (because we are on a new packet). */
        estimate_time_delay(input_items, d_n-1);
        estimate_time_delay2(input_items, d_n-1);

        estimate_channel_state(input_items, d_n-1, 0);
        interpolate_channel_state();
        for(int i=0; i<d_n ; ++i)
        {
            add_item_tag(0, nitems_written(0) + i, d_csi_key, generate_csi_pmt());
        }
        d_start_new_packet = false;
        d_correlation_offset = 1;
      } else
      {
        //std::cout << "old packet" << std::endl;
        // We are in the middle of a packet. Go on with channel estimation until the end of the packet.
        if(start_tags.size() > 0)
        {
            d_start_new_packet = true;
            items_to_process = start_tags[0].offset - nitems_read(0);
        } else
        {
            items_to_process = noutput_items;
        }
        // Iterate over OFDM symbols.
        for (int s = 0; s < items_to_process; ++s) {
          // Estimate the complex channel coefficient of all pilot carriers (of all MIMO branches).
          estimate_channel_state(input_items, s, (d_correlation_offset++)%d_n); 
          /* We have estimated the CSI for the pilot carriers.
           * Now, lets interpolate over all remaining OFDM sub-carriers. */
          interpolate_channel_state();
          /* Now we have individual CSI for each sub-carrier of each MIMO-branch.
           * Add tag with this CSI to the output vector.
           * All CSI for one time step is stored in one 3-dim vector which is tagged
           * to the output vector of the first MIMO branch. */
          add_item_tag(0, nitems_written(0) + s, d_csi_key, generate_csi_pmt());
          }
      }
            //std::cout << "items to process" << items_to_process << std::endl;

      /* Copy occupied OFDM sub-carriers to output buffer.
       * (Neither the zero-carriers nor the pilot carriers) */
      extract_payload_carriers(input_items, output_items, d_n-1, items_to_process);

      // Propagate all other tags (except the CSI tags which were changed) manually.
      std::vector <gr::tag_t> propagate_tags;
      get_tags_in_window(propagate_tags, 0, 0, items_to_process);
      for (int l = 0; l < propagate_tags.size(); ++l) {
          add_item_tag(0, propagate_tags[l].offset, propagate_tags[l].key, propagate_tags[l].value);
       //   std::cout << "propagate, tags: " << propagate_tags[l].offset << std::endl;
      }
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (items_to_process);
      // Tell runtime system how many output items we produced.
      //std::cout << "chanest written and consumed "<<items_to_process<<std::endl;
      return items_to_process;
    }

  } /* namespace digital */
} /* namespace gr */


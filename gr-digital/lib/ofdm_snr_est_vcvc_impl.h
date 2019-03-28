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

#ifndef INCLUDED_DIGITAL_OFDM_SNR_EST_VCVC_IMPL_H
#define INCLUDED_DIGITAL_OFDM_SNR_EST_VCVC_IMPL_H

#include <gnuradio/digital/ofdm_snr_est_vcvc.h>

namespace gr {
  namespace digital {

    class ofdm_snr_est_vcvc_impl : public ofdm_snr_est_vcvc
    {
     private:
      uint16_t d_num_inputs;
      uint32_t d_fft_len;
      /*! OFDM sub-carriers which are occupied with data.
       * (All not-zero and not-pilot carriers)*/
      std::vector<int> d_occupied_carriers;
      std::vector<int> d_zero_carriers;
      uint16_t d_update_time;
      float d_averaging_length;
      pmt::pmt_t d_snr_key; //!< Key for the SNR stream tags.
      std::vector<float> d_snr; /*!< Vector that stores the current SNR values of the N receivers.*/
      float *d_mag_squared;

     public:
      ofdm_snr_est_vcvc_impl(uint16_t num_inputs, uint32_t fft_len,
                             std::vector<int> occupied_carriers,
                             std::vector<int> zero_carriers,
                             const std::string &snr_key,
                             uint16_t update_time, float averaging_length);
      ~ofdm_snr_est_vcvc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_SNR_EST_VCVC_IMPL_H */


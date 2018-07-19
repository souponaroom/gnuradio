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

#ifndef INCLUDED_DIGITAL_MIMO_OFDM_CHANNEL_ESTIMATOR_VCVC_IMPL_H
#define INCLUDED_DIGITAL_MIMO_OFDM_CHANNEL_ESTIMATOR_VCVC_IMPL_H

#include <gnuradio/digital/mimo_ofdm_channel_estimator_vcvc.h>

namespace gr {
  namespace digital {

    class mimo_ofdm_channel_estimator_vcvc_impl : public mimo_ofdm_channel_estimator_vcvc
    {
     private:
      uint16_t d_n;
      uint32_t d_fft_len;
      std::vector<std::vector<gr_complex> > d_pilot_symbols;
      std::vector<int> d_pilot_carriers;
      std::vector<std::vector<std::vector<gr_complex> > > d_channel_state;

      gr_complex correlate_pilots(const gr_complex* in, std::vector<gr_complex> pilot, uint32_t distance);
      void interpolate_channel_state();

     public:
      mimo_ofdm_channel_estimator_vcvc_impl(uint16_t n,
                                            uint32_t fft_len,
                                            std::vector<std::vector<gr_complex> > pilot_symbols,
                                            std::vector<int> pilot_carriers);
      ~mimo_ofdm_channel_estimator_vcvc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_MIMO_OFDM_CHANNEL_ESTIMATOR_VCVC_IMPL_H */


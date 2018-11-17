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
      uint16_t d_n; /*!< Number of receiving and transmitting antennas. */
      uint32_t d_fft_len; /*!< FFT length = Number of OFDM sub-carriers. */
      /*!< 2-dim vector of the dimensions NxN, containing the training pilot symbols
       * for each MIMO channel. */
      std::vector<std::vector<gr_complex> > d_pilot_symbols;
      std::vector<int> d_pilot_carriers; /*!< OFDM sub-carriers, where the pilot symbols are located. */
      std::vector<int> d_occupied_carriers;
      uint32_t d_output_vlen;

      /*!< 3-dimensional vector, storing the MIMO CSI from the current OFDM symbol for
       * each sub-carrier.
       * The dimensions are fft_len, N, N.
       */
      std::vector<std::vector<std::vector<gr_complex> > > d_channel_state;

      void extract_payload_carriers(gr_vector_const_void_star &input_items,
                                    gr_vector_void_star &output_items,
                                    uint32_t length);

      void estimate_channel_state(gr_vector_const_void_star &input_items,
                                  uint32_t reading_offset);

      pmt::pmt_t generate_csi_pmt();

      /*!< Correlates the received pilot symbols with the actual ones.
       *
       * @param in Pointer to the first one of the received pilot symbols.
       * @param pilot Vector with the transmitted pilot symbols.
       * @param distance Distance between 2 pilot symbols of one sequence. In most cases, the
       * distance is equal to fft_len.
       * @param pilot_offset Shift of the pilot training sequence, to properly correlate with the first received symbol.
       * @return
       */
      gr_complex correlate_pilots(const gr_complex* in, std::vector<gr_complex> pilot, uint32_t distance, uint16_t pilot_offset);
      /*!< Linear interpolation over OFDM sub-carriers.
           The carriers beyond the edge pilot carriers hold the value of the nearest sampling points.
       */
      void interpolate_channel_state();

     public:
      mimo_ofdm_channel_estimator_vcvc_impl(uint16_t n,
                                            uint32_t fft_len,
                                            std::vector<std::vector<gr_complex> > pilot_symbols,
                                            std::vector<int> pilot_carriers,
                                            std::vector<int> occupied_carriers);
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


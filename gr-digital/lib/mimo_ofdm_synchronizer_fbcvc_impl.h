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

#ifndef INCLUDED_DIGITAL_MIMO_OFDM_SYNCHRONIZER_FBCVC_IMPL_H
#define INCLUDED_DIGITAL_MIMO_OFDM_SYNCHRONIZER_FBCVC_IMPL_H

#include <gnuradio/digital/mimo_ofdm_synchronizer_fbcvc.h>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace digital {

    class mimo_ofdm_synchronizer_fbcvc_impl : public mimo_ofdm_synchronizer_fbcvc
    {
     private:
      uint16_t d_n; /*!< Number of receiving antennas N. */
      uint32_t d_fft_len; /*!< FFT length. */
      uint32_t d_cp_len; /*!< Cyclic prefix length. */
      uint16_t d_symbol_len; /*!< FFT length + cyclic prefix length */
      /*! State variable for the synchronizer state machine. True means, that
       * we found the beginning of a frame and are currently processing it.*/
      bool d_on_frame;
      /*! Indicates first data symbol of frame (after sync symbols).
       * This is where we set a 'start' tag. */
      bool d_first_data_symbol;
      float d_phase; /*!< Phase which rotates to correct fine frequency offset. */
      int d_carrier_freq_offset; /*!< Estimated carrier frequency offset. */

      /*! The index of the first carrier with data.
       *  (index 0 is not DC here, but the lowest frequency) */
      int d_first_active_carrier;
      //! The index of the last carrier with data.
      int d_last_active_carrier;
      //! Maximum negative carrier offset. (usually a negative value!)
      int d_max_neg_carr_offset;
      //! Maximum positive carrier offset. (usually a positive value)
      int d_max_pos_carr_offset;

      fft::fft_complex *d_fft; //!< Instance of FFT class.
      std::vector<gr_complex> d_rec_sync_symbol1;
      std::vector<gr_complex> d_rec_sync_symbol2;
      std::vector<gr_complex> d_corr_v;
      pmt::pmt_t d_start_key;

      /*! \brief Rotates the phase of a complex pointer with specified params.
       * Used to correct a fractional frequency offset in OFDM.
       *
       * @param fine_freq_off Buffer with fine frequency offset for each time sample.
       * @param rotation_length Number of samples to rotate over.
       */
      void rotate_phase(const float *fine_freq_off, uint16_t rotation_length);

      /*! \brief Calculate the coarse frequency offset in number of carriers.
       * @param sync_sym1 Received synchronization symbol 1 with
       * corrected fractional frequency offset in time domain.
       * @param sync_sym2 Received synchronization symbol 2 with
       * corrected fractional frequency offset in time domain.
       * @return Carrier frequency offset. (even integer; can be negative!)
       */
      int get_carr_offset(const gr_complex *sync_sym1,
                          const gr_complex *sync_sym2);

     public:
      mimo_ofdm_synchronizer_fbcvc_impl(uint16_t n,
                                        uint32_t fft_len,
                                        uint32_t cp_len,
                                        const std::vector<gr_complex> &sync_symbol1,
                                        const std::vector<gr_complex> &sync_symbol2,
                                        const std::string &start_key);
      ~mimo_ofdm_synchronizer_fbcvc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_MIMO_OFDM_SYNCHRONIZER_FBCVC_IMPL_H */


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

#ifndef INCLUDED_DIGITAL_OFDM_CORRECT_CARRIER_FREQ_OFFSET_VCVC_IMPL_H
#define INCLUDED_DIGITAL_OFDM_CORRECT_CARRIER_FREQ_OFFSET_VCVC_IMPL_H

#include <gnuradio/digital/ofdm_correct_carrier_freq_offset_vcvc.h>

namespace gr {
  namespace digital {

    class ofdm_correct_carrier_freq_offset_vcvc_impl : public ofdm_correct_carrier_freq_offset_vcvc
    {
     private:
      uint32_t d_fft_len; /*!< Length of the FFT vectors and vlen of the items of this block. */
      uint32_t d_cp_len; /*!< Cyclic prefix length. (Required to correct the frequency offset over time. */
      std::string d_carrier_freq_offset_key; /*!< Key of the required tags which include the carrier frequency offset. */
      pmt::pmt_t d_key; /*!< PMT stores the key of the CSI tag. */
      int d_carrier_offset; /*!< Local var which remembers the current carrier frequency offset. */

      /*! \brief Writes the carrier frequency corrected fft vectors to the output.
       *
       * @param in Pointer to the input buffer.
       * @param out Poiner to the output buffer.
       * @param length Number of FFT vectors to process.
       */
      void correct_offset(const gr_complex *in, gr_complex *out, uint32_t length);

     public:
      ofdm_correct_carrier_freq_offset_vcvc_impl(uint32_t fft_len, uint32_t cp_len, std::string carrier_freq_offset_key);
      ~ofdm_correct_carrier_freq_offset_vcvc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_CORRECT_CARRIER_FREQ_OFFSET_VCVC_IMPL_H */


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


#ifndef INCLUDED_DIGITAL_MIMO_OFDM_SYNCHRONIZER_FBCVC_H
#define INCLUDED_DIGITAL_MIMO_OFDM_SYNCHRONIZER_FBCVC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace digital {

    /*! \brief MIMO-OFDM synchronization for frame, timing, fractional and integer frequency carrier offset.
     * Input ports:
     * - Trigger channel indicates the start of a frame.
     * - Frequency channel indicates fractional frequency carrier offset.
     * - Reference channel for integer frequency carrier offset estimation.
     *   (for example the sum of the MIMO RX channels)
     *
     * The following steps are applied:
     * - Wait for a trigger indicating the start of an OFDM frame.
     * - Read Schmidl & Cox synchronization symbols and estimate frequency carrier offset.
     *   (Estimated offset is tagged to stream)
     * - Copy (fractional carrier frequency corrected) data OFDM symbols to output.
     *   (Start of frame is tagged)
     * - Repeat with the next incoming start trigger.
     */
    class DIGITAL_API mimo_ofdm_synchronizer_fbcvc : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<mimo_ofdm_synchronizer_fbcvc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::mimo_ofdm_synchronizer_fbcvc.
       *
       * To avoid accidental use of raw pointers, digital::mimo_ofdm_synchronizer_fbcvc's
       * constructor is in a private implementation
       * class. digital::mimo_ofdm_synchronizer_fbcvc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint16_t n, uint32_t fft_len, uint32_t cp_len, 
                       const std::vector<gr_complex> &sync_symbol1, 
                       const std::vector<gr_complex> &sync_symbol2,
                       const std::string &start_key="start",
                       const std::string &carrier_freq_off_key="carrier_freq_offset");
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_MIMO_OFDM_SYNCHRONIZER_FBCVC_H */


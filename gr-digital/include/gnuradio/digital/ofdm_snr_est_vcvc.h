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


#ifndef INCLUDED_DIGITAL_OFDM_SNR_EST_VCVC_H
#define INCLUDED_DIGITAL_OFDM_SNR_EST_VCVC_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace digital {

    /*!
     * \brief <+description of block+>
     * \ingroup digital
     *
     */
    class DIGITAL_API ofdm_snr_est_vcvc : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<ofdm_snr_est_vcvc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of digital::ofdm_snr_est_vcvc.
       *
       * To avoid accidental use of raw pointers, digital::ofdm_snr_est_vcvc's
       * constructor is in a private implementation
       * class. digital::ofdm_snr_est_vcvc::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint16_t num_inputs, uint32_t fft_len, std::vector<int> occupied_carriers, std::vector<int> zero_carriers, const std::string &snr_key, uint16_t update_time=16, float averaging_length=8);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_SNR_EST_VCVC_H */


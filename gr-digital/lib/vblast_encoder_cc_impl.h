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

#ifndef INCLUDED_DIGITAL_VBLAST_ENCODER_CC_IMPL_H
#define INCLUDED_DIGITAL_VBLAST_ENCODER_CC_IMPL_H

#include <gnuradio/digital/vblast_encoder_cc.h>

namespace gr {
  namespace digital {
/*! \brief VBLAST encoder.
 * Encodes a serial input data stream into a vertical transmission vector
 * which is transmitted at once over multiple antennas.
 * The functionality equals a multiplexing of the data.
 */
    class vblast_encoder_cc_impl : public vblast_encoder_cc
    {
     private:
      /*! Number of output ports on which the data is divided.
       * This equals the number of your transmit antennas.*/
      uint16_t d_num_outputs; /*!< Number of transmit ports. */
      const pmt::pmt_t d_packet_len_key; /*!< PMT stores the key of the CSI tag. */

     public:
      vblast_encoder_cc_impl(uint16_t num_outputs, const std::string &packet_len_tag_key);
      ~vblast_encoder_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_VBLAST_ENCODER_CC_IMPL_H */


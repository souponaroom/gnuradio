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

#ifndef INCLUDED_DIGITAL_ALAMOUTI_DECODER_CC_IMPL_H
#define INCLUDED_DIGITAL_ALAMOUTI_DECODER_CC_IMPL_H

#include <gnuradio/digital/alamouti_decoder_cc.h>
#include <vector>

namespace gr {
  namespace digital {

    class alamouti_decoder_cc_impl : public alamouti_decoder_cc
    {
     private:
      std::vector <gr::tag_t> tags; /*!< Vector that stores the tags in input buffer. */
      static const std::string s; /*!< String that matches the key of the CSI tags. */
      static const pmt::pmt_t d_key; /*!< PMT stores the key of the CSI tag. */
      std::vector <gr_complex> d_csi;
      /*!< Array of length 2 which stores the current channel
       * state information (CSI). The array is being updated which each
       * received tag of the key='csi'.
       */
      void decode_symbol(const gr_complex* in, gr_complex* out, uint32_t length);

     public:
      alamouti_decoder_cc_impl();
      ~alamouti_decoder_cc_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_ALAMOUTI_DECODER_CC_IMPL_H */


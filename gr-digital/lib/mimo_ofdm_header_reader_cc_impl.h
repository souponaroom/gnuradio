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

#ifndef INCLUDED_DIGITAL_MIMO_OFDM_HEADER_READER_CC_IMPL_H
#define INCLUDED_DIGITAL_MIMO_OFDM_HEADER_READER_CC_IMPL_H

#include <gnuradio/digital/mimo_ofdm_header_reader_cc.h>

namespace gr {
  namespace digital {

    class mimo_ofdm_header_reader_cc_impl : public mimo_ofdm_header_reader_cc
    {
     private:
      constellation_sptr d_constellation;
      unsigned int d_dim;
      packet_header_default::sptr d_header_formatter;
      static const pmt::pmt_t d_key; /*!< PMT stores the key of the CSI tag. */
      //uint32_t d_symbol_counter;
      uint32_t d_header_length;
      uint32_t d_packet_length;
      uint32_t d_frame_length;
      uint32_t d_packet_num;
      bool d_on_packet;
      unsigned char* d_header_data;

      pmt::pmt_t d_len_tag_key;
      pmt::pmt_t d_frame_len_tag_key;
      pmt::pmt_t d_num_tag_key;

      std::vector<tag_t> d_header_tags;

      void demod_header(const gr_complex *src, unsigned char *dest);
      bool parse_header();
      void add_tags(uint32_t offset);

     public:
      mimo_ofdm_header_reader_cc_impl(constellation_sptr constellation,
                                      const gr::digital::packet_header_default::sptr &header_formatter);
      ~mimo_ofdm_header_reader_cc_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_MIMO_OFDM_HEADER_READER_CC_IMPL_H */


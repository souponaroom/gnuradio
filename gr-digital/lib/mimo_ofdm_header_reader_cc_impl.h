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
    /*! \brief Parses header and adds header tags to extracted payload stream.
     *
     * -Input: Synchronized and equalized sample stream with 'start' tags on the beginning of each frame.
     * -Read and demodulate header. Parse and validate header info.
     * -Add header info with stream tags to the beginning of each packet.
     * -Output: Payload stream (not yet demodulated) with header tags at the beginning of each packet.
     */

    class mimo_ofdm_header_reader_cc_impl : public mimo_ofdm_header_reader_cc
    {
     private:
      //! Constellation object, defining the applied digital constellation.
      constellation_sptr d_constellation;
      //! Dimensionality of the applied digital constellation.
      unsigned int d_dim;
      //! Header formatter object.
      packet_header_default::sptr d_header_formatter;
      pmt::pmt_t d_start_key; //!< Key of the start stream tags.
      //! Length of header in number of complex symbols.
      uint32_t d_header_length;
      //! Payload length of packet in number of complex symbols (without header length).
      uint32_t d_packet_length;
      //! Length of a frame in number of complex symbols (1 frame can contain multiple packets).
      uint32_t d_frame_length;
      //! Identification number of packet (usually counting up from 0).
      uint32_t d_packet_num;
      //! Indicates whether we are currently processing a packet.
      bool d_on_packet;
      //! Buffer in which we write the demodulated header.
      unsigned char* d_header_data;

      uint32_t d_remaining_packet_len;

      /*! PMT stores the key of the packet (payload) length tag. */
      pmt::pmt_t d_len_tag_key;
      /*! PMT stores the key of the frame length tag. */
      pmt::pmt_t d_frame_len_tag_key;
      /*! PMT stores the key of the identification number tag. */
      pmt::pmt_t d_num_tag_key;

      //! Vector which stores all extracted tags from the current packet header.
      std::vector<tag_t> d_header_tags;

      uint32_t locate_next_tag(uint32_t offset, uint32_t buffer_length);

      /*! \brief Demodulation of the header symbols.
       * Decides for the constellation point with the minimal distance
       * to the sample, using the Euclidean distance as a metric.
       *
       * @param src Input buffer with received header samples.
       * @param dest Output buffer where the demodulated symbols are written to.
       */
      void demod_header(const gr_complex *src, unsigned char *dest);

      /*! \brief Check header for validity and parse header info.
       *
       * @return True if header is valid, false if not.
       */
      bool parse_header();
      /*! \brief Add stream tags with header info to payload stream.
       *
       * @param offset Position of the added tags in the output buffer.
       */
      void add_tags(uint32_t offset);

     public:
      mimo_ofdm_header_reader_cc_impl(constellation_sptr constellation,
                                      const gr::digital::packet_header_default::sptr &header_formatter,
                                      const std::string &start_key, const std::string &len_tag_key,
                                      const std::string &frame_len_tag_key, const std::string &num_tag_key);
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

